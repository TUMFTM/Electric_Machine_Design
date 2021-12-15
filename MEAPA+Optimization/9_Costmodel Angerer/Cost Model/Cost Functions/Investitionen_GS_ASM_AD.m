%Definition of the investment costs for an ASM for large-scale production.
%Detailed source information for each individual machine
%Data in €

%% Machine costs (unit prices)

Druckgussmaschine=1500000;      
Entgratpresse=150000;           
BAZ=400000;
Reinigungsanlage=133000;
Stanzanlage=1500000;
Paketieranlage=3000000;
Isolieranlage=200000;
Wickelautomat=2725000;
Bandagieranlage=200000;
Pruefanlage_ST=300000;
Impraegrnieranlage=800000;
Umformpresse=700000;
Reinigunganlage=150000;
Verkettung=1643400;
Induktionsanlage=250000;
Unwucht_Testanlage=470000;
Magnetbestueckung=860000;
Montagestation=150000;
Statormontage=900000;
Rotormontage=150000;
Anbauteilmontage=100000;
Endpruefung=1500000;

t_EM=190;

%% Housing production

t_DG=100; %cycle time die casting
t_EG=60; %cycle time trimming press
t_BAZ=180; %cycle time machining
t_RE=250; %cycle time cleaning
t_VK=180; %cycle time linking

Anzahl_Maschinen_DG_GH=ceil(t_DG/t_EM);
Anzahl_Maschinen_EG_GH=ceil(t_EG/t_EM);   
Anzahl_Maschinen_BAZ_GH=ceil(t_BAZ/t_EM); 
Anzahl_Maschinen_RE_GH=ceil(t_RE/t_EM); 
Anzahl_Maschinen_VK_GH=ceil(t_VK/t_EM);

Druckguss_GH=Druckgussmaschine*Anzahl_Maschinen_DG_GH;
Entgraten_GH=Entgratpresse*Anzahl_Maschinen_EG_GH;
BAZ_GH=BAZ*Anzahl_Maschinen_BAZ_GH;
Reinigung_GH=Reinigungsanlage*Anzahl_Maschinen_RE_GH;
Verkettung_GH=Verkettung*Anzahl_Maschinen_VK_GH;

Invest_GH=Druckguss_GH+Entgraten_GH+BAZ_GH+Reinigung_GH+Verkettung_GH;

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
t_WI=260; %cycle time winding
t_BAN=320; %cycle time machining
t_PR=60; %Takt time testing
t_IM=85; %cycle time impregnation
t_VK=180; %cycle time linking

Anzahl_Maschinen_NI_ST=ceil(t_NI/t_EM);
Anzahl_Maschinen_WI_ST=ceil(t_WI/t_EM);   
Anzahl_Maschinen_BAN_ST=ceil(t_BAN/t_EM); 
Anzahl_Maschinen_PR_ST=ceil(t_PR/t_EM); 
Anzahl_Maschinen_IM_ST=ceil(t_IM/t_EM);
Anzahl_Maschinen_VK_ST=ceil(t_VK/t_EM);

Nutisolierung_ST=Isolieranlage*Anzahl_Maschinen_NI_ST;
Spulenwicklung_ST=Wickelautomat*Anzahl_Maschinen_WI_ST;
Bandagierung_ST=Bandagieranlage*Anzahl_Maschinen_BAN_ST;
Test_ST=Pruefanlage_ST*Anzahl_Maschinen_PR_ST;
Impraegnierung_ST=Impraegrnieranlage*Anzahl_Maschinen_IM_ST;
Verkettung_ST=Verkettung*Anzahl_Maschinen_VK_ST;

Invest_ST=Nutisolierung_ST+Spulenwicklung_ST+Bandagierung_ST+Test_ST+Impraegnierung_ST+Verkettung_ST;

%% Shaft production

t_FG=60; %cycle time forming
t_RH=40; %cycle time surface hardening
t_BAZ=320; %cycle time machining
t_RE=180; %Cleaning cycle time
t_VK=200; %cycle time interlinking

Anzahl_Maschinen_FG_W=ceil(t_FG/t_EM);
Anzahl_Maschinen_RH_W=ceil(t_RH/t_EM);   
Anzahl_Maschinen_BAZ_W=ceil(t_BAZ/t_EM); 
Anzahl_Maschinen_RE_W=ceil(t_RE/t_EM); 
Anzahl_Maschinen_VK_W=ceil(t_VK/t_EM);

Formgebung_W=Umformpresse*Anzahl_Maschinen_FG_W;
Haerten_W=Induktionsanlage*Anzahl_Maschinen_RH_W;
BAZ_W=BAZ*Anzahl_Maschinen_BAZ_W;
Reinigung_W=Reinigunganlage*Anzahl_Maschinen_RE_W;
Verkettung_W=Verkettung*Anzahl_Maschinen_VK_W;

Invest_W=Formgebung_W+Haerten_W+BAZ_W+Reinigung_W+Verkettung_W;

%% Rotor production

t_DG=110; %Takt time magnet assembly
t_WM=300; %cycle time shaft assembly
t_UW=120; %cycle time unbalance test
t_VK=200; %cycle time interlinking

Anzahl_Maschinen_MG_RT=ceil(t_DG/t_EM);
Anzahl_Maschinen_WM_RT=ceil(t_WM/t_EM);   
Anzahl_Maschinen_UW_RT=ceil(t_UW/t_EM); 
Anzahl_Maschinen_VK_RT=ceil(t_VK/t_EM);

Druckguss_RT=Druckgussmaschine*Anzahl_Maschinen_MG_RT;
Wellenmontage_RT=Montagestation*Anzahl_Maschinen_WM_RT;
Unwucht_RT=Unwucht_Testanlage*Anzahl_Maschinen_UW_RT;
Verkettung_RT=Verkettung*Anzahl_Maschinen_VK_RT;

Invest_RT=Druckguss_RT+Wellenmontage_RT+Unwucht_RT+Verkettung_RT;

%% Final assembly

t_SM=200; %Takt time stator assembly
t_RM=500; %Takt time rotor assembly
t_AM=850; %Takt time mounting part assembly
t_ET=160; %cycle time final inspection

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
    Pruefanlage_ST Reinigunganlage Reinigungsanlage Rotormontage...
    Stanzanlage Statormontage t_AM t_BAN t_BAZ t_DG t_EG t_EM t_ET t_FG...
    t_IM t_MG t_NI t_PA t_PR t_RE t_RH t_RM t_SM t_ST t_UW t_VK t_WI...
    t_WM Druckguss_RT Umformpresse Unwucht_Testanlage Verkettung Wickelautomat