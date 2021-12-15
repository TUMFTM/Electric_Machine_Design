%Definition of investment costs for a gearbox large series production
%Gearbox 2-stage 1-speed
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
Bearbeitungszentetrum=350000;
Impraegrnieranlage=800000;
Umformpresse=700000;
Reinigunganlage=150000;
Waelzfraesmaschine=500000;
Raeummaschine=581000;
Schleifmaschine=600000;
Verkettung=1000000;
Induktionsanlage=250000;
Unwucht_Testanlage=470000;
Magnetbestueckung=860000;
Montagestation=150000;
Statormontage=900000;
Rotormontage=150000;
Zahnradmontage=50000;
Lagermontage=50000;
Endtest_GT=400000;
Anbauteilmontage=100000;
Endpruefung=1500000;

t_EM=190; %Takt time of gearbox production

%% Housing production

t_DG=100; %cycle time die casting
t_EG=60; %cycle time trimming press
t_BAZ=180; %cycle time machining
t_RE=250; %cycle time cleaning

Anzahl_Maschinen_DG_GH=ceil(t_DG/t_EM);
Anzahl_Maschinen_EG_GH=ceil(t_EG/t_EM);   
Anzahl_Maschinen_BAZ_GH=ceil(t_BAZ/t_EM); 
Anzahl_Maschinen_RE_GH=ceil(t_RE/t_EM); 

Druckguss_GH=Druckgussmaschine*Anzahl_Maschinen_DG_GH;
Entgraten_GH=Entgratpresse*Anzahl_Maschinen_EG_GH;
BAZ_GH=BAZ*Anzahl_Maschinen_BAZ_GH;
Reinigung_GH=Reinigungsanlage*Anzahl_Maschinen_RE_GH;


Invest_GH=Druckguss_GH+Entgraten_GH+BAZ_GH+Reinigung_GH;

%% Shaft production

t_FG=60; %cycle time forming
t_RH=40; %cycle time surface hardening
t_BAZ=750; %cycle time machining
t_RE=180; %Cleaning cycle time
t_UW=100; %cycle time interlinking

Anzahl_Maschinen_FG_W=ceil(t_FG/t_EM);
Anzahl_Maschinen_RH_W=ceil(t_RH/t_EM);   
Anzahl_Maschinen_BAZ_W=ceil(t_BAZ/t_EM); 
Anzahl_Maschinen_RE_W=ceil(t_RE/t_EM); 
Anzahl_Maschinen_UW_W=ceil(t_UW/t_EM);

Formgebung_W=Umformpresse*Anzahl_Maschinen_FG_W;
Haerten_W=Induktionsanlage*Anzahl_Maschinen_RH_W;
BAZ_W=Bearbeitungszentetrum*Anzahl_Maschinen_BAZ_W;
Reinigung_W=Reinigunganlage*Anzahl_Maschinen_RE_W;
Unwucht_W=Unwucht_Testanlage*Anzahl_Maschinen_UW_W;

Invest_W=Formgebung_W+Haerten_W+BAZ_W+Reinigung_W+Unwucht_W;

%% Gear manufacturing

t_FG=60; %Takt time shaping
t_WF=1880; %cycle time gear hobbing
t_RM=180; %Broaching cycle time
t_RE=320; %cycle time cleaning
t_HA=100; %cycle time induction hardening
t_SL=750; %cycle time grinding

Anzahl_Maschinen_FG_ZR=ceil(t_FG/t_EM);
Anzahl_Maschinen_WF_ZR=ceil(t_WF/t_EM);   
Anzahl_Maschinen_RM_ZR=ceil(t_RM/t_EM); 
Anzahl_Maschinen_RE_ZR=ceil(t_RE/t_EM); 
Anzahl_Maschinen_HA_ZR=ceil(t_HA/t_EM);
Anzahl_Maschinen_SL_ZR=ceil(t_SL/t_EM);

Formgebung_ZR=Umformpresse*Anzahl_Maschinen_FG_ZR;
Waelzfraesen_ZR=Waelzfraesmaschine*Anzahl_Maschinen_WF_ZR;
Rauemen_ZR=Raeummaschine*Anzahl_Maschinen_RM_ZR;
Reinigung_ZR=Reinigungsanlage*Anzahl_Maschinen_RE_ZR;
Haerten_ZR=Induktionsanlage*Anzahl_Maschinen_HA_ZR;
Schleifen_ZR=Schleifmaschine*Anzahl_Maschinen_SL_ZR;

Invest_ZR=Formgebung_ZR+Waelzfraesen_ZR+Rauemen_ZR+Reinigung_ZR+Haerten_ZR+Schleifen_ZR;

%% Final assembly

t_ZR=160; %Takt time gear assembly
t_LA=160; %cycle time bearing assembly
t_AM=420; %Takt time mounting part assembly
t_ET=160; %cycle time final inspection

Anzahl_Maschinen_ZR_EM=ceil(t_ZR/t_EM);
Anzahl_Maschinen_LA_EM=ceil(t_LA/t_EM);   
Anzahl_Maschinen_AM_EM=ceil(t_AM/t_EM); 
Anzahl_Maschinen_ET_EM=ceil(t_ET/t_EM);

ZRmontage_EM=Zahnradmontage*Anzahl_Maschinen_ZR_EM;
Lagermontage_EM=Lagermontage*Anzahl_Maschinen_LA_EM;
Anbauteilmontage_EM=Anbauteilmontage*Anzahl_Maschinen_AM_EM;
Endtest_EM=Endtest_GT*Anzahl_Maschinen_ET_EM;

Invest_EM=ZRmontage_EM+Lagermontage_EM+Anbauteilmontage_EM+Endtest_EM;

%% Total investment costs

Invest=Invest_GH+Invest_W+Invest_ZR+Invest_EM;

clearvars Anbauteilmontage_EM Bandagierung_ST BAZ_BP BAZ_GH BAZ_W...
    Druckguss_GH Druckguss_RT Endtest_EM Entgraten_GH Formgebung_W...
    Haerten_W Impraegnierung_ST Invest_BP Invest_EM Invest_GH Invest_RT...
    Invest_ST Invest_W Nutisolierung_ST Paketieranlage_BP Reinigung_GH...
    Reinigung_W Rotormontage_EM Spulenwicklung_ST Stanzen_BP...
    Statormontage_EM Test_ST Unwucht_RT Verkettung_GH Verkettung_RT...
    Verkettung_ST Verkettung_W Wellenmontage_RT Formgebung_ZR Haerten_ZR...
    Invest_ZR Lagermontage_EM Rauemen_ZR Reinigung_ZR Schleifen_ZR Unwucht_W...
    Waelzfraesen_ZR WZ_ZR ZRmontage_EM...
    Anbauteilmontage_EM Bandagierung_ST BAZ_BP BAZ_GH BAZ_W...
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
    t_WM Druckguss_RT Umformpresse Unwucht_Testanlage Verkettung Wickelautomat...
    Anzahl_Maschinen_FG_ZR Anzahl_Maschinen_HA_ZR Anzahl_Maschinen_LA_EM...
    Anzahl_Maschinen_RE_ZR Anzahl_Maschinen_RM_ZR Anzahl_Maschinen_SL_ZR...
    Anzahl_Maschinen_UW_W Anzahl_Maschinen_WF_ZR Anzahl_Maschinen_ZR_EM...
    Bearbeitungszentetrum Endtest_GT Lagermontage Raeummaschine...
    Schleifmaschine t_HA t_LA t_SL t_WF t_ZR Waelzfraesmaschine Zahnradmontage