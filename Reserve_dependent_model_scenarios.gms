*** MODEL FOR ONSHORE AND OFFSHORE MINERAL EXTRACTION ***



*** KEY INSTRUCTIONS ***

*** STEP 1: Check model settings
*** STEP 2: Check that time horizon and step-size align
*** STEP 3: Activate desired functions (depending on whether interested in monopoly or duopoly solution) (see in-code description)
***             - If interested in the duopoly solution:
***                 (1) Max obj. of sector A subject to zero-decision vector for sector B. Save the results in the putfile (results.dat) to the csv file (Input_data_for_finding_cournot-nash-equilibrium.csv). 
***                 (2) Max obj. of sector B subject to decision vector for sector A from solution in (1).
***                 (3) Repeat the process in (2) while alternating between the sectors until the decision vectors stop changing from one solution to the next
***                 (4) Do robustness tests (e.g. alternate which sector you start with, and start with decision vector different from a zero-vector. Check that the results converge towards the same solution.)
*** STEP 4: Store the final results from results.dat in the DATA_TEMPLATE file.
*** STEP 5: Open Code_for_plotting_results.r in R, and import data from DATA_TEMPLATE.xlsx, and run the R-code to plot the results. 



*** SOME OTHER NOTES ***

* CHOOSE TIME HORIZON CAREFULLY (LONGER TIME HORIZON AND POSITIVE DISCOUNT RATE WILL REQUIRE HIGH SCALING OF OBJ / LOW TOLERANCE FOR CONVERGENCE)
* USE APPROPRIATE TIME STEP (AS HIGH AS POSSIBLE WHILE STILL AVOIDING NUMERICAL INTEGRATION ERROR)
* USE PROPER SCALING OF OBJECTIVE (AS LOW AS POSSIBLE WHILE STILL ACHIEVING ROBUST RESULTS)
* ALTERNATIVE TO THE ABOVE: PROPER ADJUSTMENT OF TOLERANCE FOR CONVERGENCE (AS HIGH AS POSSIBLE WHILE STILL ACHIEVING ROBUST RESULTS)
* SMALL STEP-SIZE WILL INCREASE SOLVING TIME, BUT MAY BE NEEDED TO AVOID SIGNIFICANT NUMERICAL INTEGRATION ERROR
* UPWARD SCALING OF OBJECTIVE / DOWNWARD SCALING OF TOLERANCE WILL INCREASE SOLVING TIME, BUT MAY BE NEEDED, ESPECIALLY FOR ACCURATE RESULTS FAR OUT IN THE TIME HORIZON
* APPROPRIATE TIME STEP AND SCALING DEPENDS ON SCENARIO AND NUMERICAL SPECIFICATIONS...
* GOOD INITIAL GUESSES FOR CONTROL VARIABLES SAVES TIME ON THE SOLVING PROCESS


Option DECIMALS=8
Option solver=knitro
Option iterlim = 20000


SET
T TIME /T-1*T-2000/     
I RESOURCE /LAND, HAV/
J LEVELVAR /U1, U2, INV1, INV2, K1, K2, X1, X2, NPVTOT, NPVLAND, NPVOCEAN/;

*NOTE: MAKE SURE FINAL TIME = FINAL YEAR/TIME STEP (CONVERT BACK WHEN PLOTTING RESULTS

SCALAR
PMAX        PRICE WHEN PRODUCTION IS ZERO  /0.006/
DISK        DISCOUNT RATE                  /0.05/
PC          PRICE CURVATURE PARAMETER      /0.0001/
H           STEP-SIZE                      /0.1/
SF          OBJ VAL SCALING FACTOR         /1e2/;


PARAMETERS
A(I)        TOTAL FACTOR PRODUCTIVITY
/LAND 0.0003
HAV   0.0001/

D(I)        DEPRECIATION RATE
/LAND 0.1
HAV   0.1/


C(I)        Cost Parameter
/LAND 0.3
HAV   0.3/

F(I)        Cost Parameter
/LAND 0.5
HAV   0.5/


E(I)        Cost Parameter
/LAND 1.1
HAV   1.1/


Table g(T,J) RESULTS
$ondelim
$include data.csv
$offdelim;

VARIABLES
NPV         NET PRESENT VALUE
NPVT(T)     NET PRESENT VALUE TRACKER
NPVT1(T)     NET PRESENT VALUE TRACKER LAND
NPVT2(T)     NET PRESENT VALUE TRACKER OCEAN
PI(T)       PROFIT RATE
X(I,T)      RESOURCE RESERVES
K(I,T)      PRODUCTION CAPITAL
INV(I,T)    INVESTMENT
U(I,T)      Production;


POSITIVE VARIABLES
X, K, INV;

** INITIAL VALUES **
K.FX('LAND','T-1') = 40;
K.FX('HAV','T-1')  = 0;
X.FX('LAND','T-1') = 2000000;
X.FX('HAV','T-1')  = 3000000;
NPVT.FX('T-1') = 0;
NPVT1.FX('T-1') = 0;
NPVT2.FX('T-1') = 0;

** INITIAL GUESSES TO HELP SOLVER OR RUN TESTS **
*INV.L('HAV', T)=g(T,'INV2');
*U.L('HAV',T)=g(T,'U2');
*INV.L('LAND',T)=g(T,'INV1');
*U.L('LAND',T)=g(T, 'U1');

*INV.L('HAV', T)=0;
*U.L('HAV',T)=0;
*INV.L('LAND',T)=0;
*U.L('LAND',T)=0;



****** DUOPOLY CONSTRAINTS *******

**** ACTIVATE IN FIRST RUN.... 
*INV.FX('HAV', T)=0;
*U.FX('HAV',T)=0;

*INV.FX('LAND',T)=0;
*U.FX('LAND',T)=0;



**** ACTIVATE WHEN MAX LAND GIVEN WHAT OCEAN DOES (REMEMBER TO UPDATE CSV FILE)
*INV.FX('HAV', T)=g(T,'INV2');
*U.FX('HAV',T)=g(T,'U2');


**** ACTIVATE WHEN MAX OCEAN GIVEN WHAT LAND DOES (REMEMBER TO UPDATE CSV FILE)
*U.FX('LAND',T)=g(T, 'U1');
*INV.FX('LAND',T)=g(T,'INV1');



EQUATIONS
** DECLARING EQUATIONS ** 
OBJ
NPVTID
NPVTID1
NPVTID2
NR(T)
DYNMIN(I,T)
DYNCAP(I,T)
CON1(I,T)
CON2(I,T)
CON3(I,T)
CON4(I,T)
CON5(I,T);

** OBJECTIVE FUNCTION - DISCOUNTED NET PROFITS **
*OBJ..              NPV =E=  NPVT('T-100');
OBJ..               NPV =E=  SUM(T,(1/(1+DISK))**((ORD(T)-1)*H)*PI(T))*SF;


** DEFINING PROFIT RATE FUNCTION ENTERING OBJECTIVE FUNCTION ABOVE **

***ORIGINAL*** ACTIVATE WHEN MONOPOLY
NR(T)..             PI(T) =E= PMAX/(1+PC*sum(I,X(I,T)-X(I,T+1))/H)*sum(I,X(I,T)-X(I,T+1)) - sum(I,(C(I)*U(I,T))/(1+X(I,T)*A(I)))*H - (sum(I,F(I)*INV(I,T)**E(I)))*H;



***DUOPOLY***

********LAND******* ACTIVATE WHEN MAX LAND GIVEN WHAT OCEAN DOES
*NR(T)..             PI(T) =E= PMAX/(1+PC*sum(I,X(I,T)-X(I,T+1))/H)*(X('LAND',T)-X('LAND',T+1)) - C('LAND')*U('LAND',T)/(1+X('LAND',T)*A('LAND'))*H - F('LAND')*INV('LAND',T)**E('LAND')*H;

********HAV******** ACTIVATE WHEN MAX OCEAN GIVEN WHAT LAND DOES
*NR(T)..             PI(T) =E= PMAX/(1+PC*sum(I,X(I,T)-X(I,T+1))/H)*(X('HAV',T)-X('HAV',T+1)) - C('HAV')*U('HAV',T)/(1+X('HAV',T)*A('HAV'))*H - F('HAV')*INV('HAV',T)**E('HAV')*H;



** DYNAMICS - EULER APPROACH **
DYNCAP(I,T+1)..     K(I,T+1) =E= K(I,T)- K(I,T)*D(I)*H+INV(I,T)*H;
DYNMIN(I,T+1)..     X(I,T+1) =E= X(I,T) - U(I,T)*H;
NPVTID1(T+1)..      NPVT1(T+1) =E= NPVT1(T) + (1/(1+DISK))**((ORD(T)-1)*H)*(PMAX/(1+PC*sum(I,X(I,T)-X(I,T+1))/H)*(X('LAND',T)-X('LAND',T+1)) - C('LAND')*U('LAND',T)/(1+X('LAND',T)*A('LAND'))*H - F('LAND')*INV('LAND',T)**E('LAND')*H);
NPVTID2(T+1)..      NPVT2(T+1) =E= NPVT2(T) + (1/(1+DISK))**((ORD(T)-1)*H)*(PMAX/(1+PC*sum(I,X(I,T)-X(I,T+1))/H)*(X('HAV',T)-X('HAV',T+1)) - C('HAV')*U('HAV',T)/(1+X('HAV',T)*A('HAV'))*H - F('HAV')*(INV('HAV',T)**E('HAV'))*H);
NPVTID(T+1)..       NPVT(T+1) =E= NPVT1(T) + (1/(1+DISK))**((ORD(T)-1)*H)*(PMAX/(1+PC*sum(I,X(I,T)-X(I,T+1))/H)*(X('LAND',T)-X('LAND',T+1)) - C('LAND')*U('LAND',T)/(1+X('LAND',T)*A('LAND'))*H - F('LAND')*INV('LAND',T)**E('LAND')*H) + NPVT2(T) + (1/(1+DISK))**((ORD(T)-1)*H)*(PMAX/(1+PC*sum(I,X(I,T)-X(I,T+1))/H)*(X('HAV',T)-X('HAV',T+1)) - C('HAV')*U('HAV',T)/(1+X('HAV',T)*A('HAV'))*H - F('HAV')*INV('HAV',T)**E('HAV')*H);

** SPECIAL CONSTRAINTS **
CON1(I,T)..        U(I,T) =L= A(I)*K(I,T)*X(I,T);


CON2(I,T)..        U('LAND',T) =G= 0;
CON3(I,T)..        U('HAV',T) =G= 0;


CON4(I,T)..        INV('LAND',T) =G= 0;
CON5(I,T)..        INV('HAV',T) =G= 0;

MODEL
TOTAL /ALL/;


SOLVE
TOTAL USING NLP MAXIMIZING NPV;


DISPLAY X.L, K.L, INV.L, NPVT.L;


file results /M:rasmus.dat/;
results.nd=8


put results;


put /' ',',','U1',',', 'U2', ',','INV1',',', 'INV2', ',' 'K1',',','K2',',', 'X1', ',', 'X2',',', 'NPVTOT', ',', 'NPVLAND', ',', 'NPVOCEAN'/
loop((t), put t.tl,',',U.l('land',t), ',',U.l('hav',t),',',inv.l('land',t), ',',inv.l('hav',t),','k.l('land',t),',',k.l('hav',t), ',', x.l('land',t),',',x.l('hav',t), ',', NPVT.L(T), ',', NPVT1.L(T),',',NPVT2.L(T)/);

