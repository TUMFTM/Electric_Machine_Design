%Definition of the investment costs for an ASM for prototype production.
%Detailed sources for each individual machine
%Data in €

%% Machine costs (unit prices)

Druckgussmaschine=1000000;      
Entgratpresse=150000;           
BAZ=400000;
Reinigungsanlage=58000;
Stanzanlage=100000;
Paketieranlage=100000;
Runddrehen=180000;
Isolieranlage=15000;
Wickelautomat=500000;
Bandagieranlage=250000;
Pruefanlage_ST=10000;
Impraegrnieranlage=15000;
Umformpresse=700000;
Reinigunganlage=150000;
Verkettung=1000000;
Induktionsanlage=170000;
Unwucht_Testanlage=100000;
Magnetbestueckung=860000;
Montagestation=150000;
Statormontage=100000;
Rotormontage=50000;
Anbauteilmontage=10000;
Endpruefung=45000;

t_EM=190;

%% Housing production

t_DG=100; %cycle time die casting
t_BAZ=180; %cycle time machining
t_RE=250; %cycle time cleaning

Anzahl_Maschinen_DG_GH=ceil(t_DG/t_EM);
Anzahl_Maschinen_BAZ_GH=ceil(t_BAZ/t_EM); 
Anzahl_Maschinen_RE_GH=ceil(t_RE/t_EM); 

Druckguss_GH=Druckgussmaschine*Anzahl_Maschinen_DG_GH;
BAZ_GH=Runddrehen*Anzahl_Maschinen_BAZ_GH;
Reinigung_GH=Reinigungsanlage*Anzahl_Maschinen_RE_GH;

Invest_GH=Druckguss_GH+BAZ_GH+Reinigung_GH;

%% Sheet bundle creation

t_ST=60; %Punching cycle time
t_PA=180; %Packaging cycle time
t_BAZ=180; %cycle time machining

Anzahl_Maschinen_ST_BP=ceil(t_ST/t_EM);
Anzahl_Maschinen_PA_BP=ceil(t_PA/t_EM);   
Anzahl_Maschinen_BAZ_BP=ceil(t_BAZ/t_EM); 

Stanzen_BP=Stanzanlage*Anzahl_Maschinen_ST_BP;
Paketieranlage_BP=Paketieranlage*Anzahl_Maschinen_PA_BP;
BAZ_BP=BAZ*Anzahl_Maschinen_BAZ_BP;

Invest_BP=Stanzen_BP+Paketieranlage_BP+BAZ_BP;

%% Stator production

t_NI=160; %Takt time slot insulation
t_WI=170; %cycle time winding
t_BAN=320; %cycle time machining
t_PR=60; %Takt time testing
t_IM=85; %cycle time impregnation


Anzahl_Maschinen_NI_ST=ceil(t_NI/t_EM);
Anzahl_Maschinen_WI_ST=ceil(t_WI/t_EM);   
Anzahl_Maschinen_BAN_ST=ceil(t_BAN/t_EM); 
Anzahl_Maschinen_PR_ST=ceil(t_PR/t_EM); 
Anzahl_Maschinen_IM_ST=ceil(t_IM/t_EM);


Nutisolierung_ST=Isolieranlage*Anzahl_Maschinen_NI_ST;
Spulenwicklung_ST=Wickelautomat*Anzahl_Maschinen_WI_ST;
Bandagierung_ST=Bandagieranlage*Anzahl_Maschinen_BAN_ST;
Test_ST=Pruefanlage_ST*Anzahl_Maschinen_PR_ST;
Impraegnierung_ST=Impraegrnieranlage*Anzahl_Maschinen_IM_ST;


Invest_ST=Nutisolierung_ST+Spulenwicklung_ST+Bandagierung_ST+Test_ST+Impraegnierung_ST;

%% Shaft production


t_RH=40; %cycle time surface hardening
t_BAZ=180; %cycle time machining
t_RE=180; %cycle time cleaning

Anzahl_Maschinen_RH_W=ceil(t_RH/t_EM);   
Anzahl_Maschinen_BAZ_W=ceil(t_BAZ/t_EM); 
Anzahl_Maschinen_RE_W=ceil(t_RE/t_EM); 

Haerten_W=Induktionsanlage*Anzahl_Maschinen_RH_W;
BAZ_W=BAZ*Anzahl_Maschinen_BAZ_W;
Reinigung_W=Reinigunganlage*Anzahl_Maschinen_RE_W;

Invest_W=Haerten_W+BAZ_W+Reinigung_W;

%% Rotor production

t_DG=200; %Takt time magnet assembly
t_WM=300; %cycle time shaft assembly
t_UW=120; %cycle time unbalance test

Anzahl_Maschinen_MG_RT=ceil(t_DG/t_EM);
Anzahl_Maschinen_WM_RT=ceil(t_WM/t_EM);   
Anzahl_Maschinen_UW_RT=ceil(t_UW/t_EM); 

Druckguss_RT=Rotormontage*Anzahl_Maschinen_MG_RT;
Wellenmontage_RT=Montagestation*Anzahl_Maschinen_WM_RT;
Unwucht_RT=Unwucht_Testanlage*Anzahl_Maschinen_UW_RT;

Invest_RT=Druckguss_RT+Wellenmontage_RT+Unwucht_RT;

%% Final assembly

t_SM=160; %Takt time stator assembly
t_RM=300; %Takt time rotor assembly
t_AM=850; %Takt time mounting part assembly
t_ET=160; %Takt time final inspection

Anzahl_Maschinen_SM_EM=ceil(t_SM/t_EM);
Anzahl_Maschinen_RM_EM=ceil(t_RM/t_EM);   
Anzahl_Maschinen_AM_EM=ceil(t_AM/t_EM); 
Anzahl_Maschinen_ET_EM=ceil(t_ET/t_EM);

Statormontage_EM=Statormontage*Anzahl_Maschinen_SM_EM;
Rotormontage_EM=Rotormontage*Anzahl_Maschinen_RM_EM;
Anbauteilmontage_EM=Anbauteilmontage*Anzahl_Maschinen_AM_EM;
Endtest_EM=Endpruefung*Anzahl_Maschinen_ET_EM;

Invest_EM=Statormontage_EM+Rotormontage_EM+Anbauteilmontage_EM+Endtest_EM;

%% Total investment costs

Invest=Invest_GH+Invest_BP+Invest_ST+Invest_W+Invest_RT+Invest_EM;

clearvars Anbauteilmontage_EM Bandagierung_ST BAZ_BP BAZ_GH BAZ_W...
    Druckguss_GH Magnetbest_RT Endtest_EM Entgraten_GH Formgebung_W...
    Haerten_W Impraegnierung_ST Invest_BP Invest_EM Invest_GH Invest_RT...
    Invest_ST Invest_W Nutisolierung_ST Paketieranlage_BP Reinigung_GH...
    Reinigung_W Rotormontage_EM Spulenwicklung_ST Stanzen_BP...
    Statormontage_EM Test_ST Unwucht_RT Verkettung_GH Verkettung_RT...
    Verkettung_ST Verkettung_W Wellenmontage_RT Anbauteilmontage Anzahl_Maschinen_AM_EM Anzahl_Maschinen_BAN_ST...
    Anzahl_Maschinen_BAZ_BP Anzahl_Maschinen_BAZ_GH Anzahl_Maschinen_BAZ_W...
    Anzahl_Maschinen_DG_GH Anzahl_Maschinen_EG_GH Anzahl_Maschinen_ET_EM...
    Anzahl_Maschinen_FG_W Anzahl_Maschinen_IM_ST Anzahl_Maschinen_MG_RT...
    Anzahl_Maschinen_NI_ST Anzahl_Maschinen_PA_BP Anzahl_Maschinen_PR_ST...
    Anzahl_Maschinen_RE_GH Anzahl_Maschinen_RE_W Anzahl_Maschinen_RH_W...
    Anzahl_Maschinen_RM_EM Anzahl_Maschinen_SM_EM Anzahl_Maschinen_ST_BP...
    Anzahl_Maschinen_UW_RT Anzahl_Maschinen_VK_GH Anzahl_Maschinen_VK_RT...
    Anzahl_Maschinen_VK_ST Anzahl_Maschinen_VK_W Anzahl_Maschinen_WI_ST...
    Anzahl_Maschinen_WM_RT Bandagieranlage BAZ Druckgussmaschine...
    Endpruefung Entgratpresse Impraegrnieranlage Induktionsanlage...
    Isolieranlage Magnetbestueckung Montagestation Paketieranlage ...
    Pruefanlage_ST Runddrehen Druckguss_RT Reinigunganlage Reinigungsanlage Rotormontage...
    Stanzanlage Statormontage t_AM t_BAN t_BAZ t_DG t_EG t_EM t_ET t_FG...
    t_IM t_MG t_NI t_PA t_PR t_RE t_RH t_RM t_SM t_ST t_UW t_VK t_WI...
    t_WM Umformpresse Unwucht_Testanlage Verkettung Wickelautomat