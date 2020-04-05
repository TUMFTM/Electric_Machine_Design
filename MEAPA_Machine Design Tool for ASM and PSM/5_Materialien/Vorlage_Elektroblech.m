% -------------------------------------------------------------------------
% TU Muenchen - Lehrstuhl fuer Fahrzeugtechnik (FTM)
% -------------------------------------------------------------------------
% Modell fuer den Entwurf und die Analyse einer PMSM oder ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor: Svenja Kalt (kalt@ftm.mw.tum.de)
%        Jonathan Erhard
% -------------------------------------------------------------------------

% Hinweis: Um neue Dateien fuer Elektrobleche zu erzeugen muss ein struct 
%          mit den Eintraegen
%           - Bezeichnung
%           - rho_Fe
%           - H
%           - B
%           - mu_r
%           - p_vFe_B_vec
%           - p_vFe_f_vec
%           - p_vFe(Spalten: p_vFe_f_vec, Zeilen: p_vFe_B_vec)
%          angelegt werden (siehe VACOFLUX 48). Anschliessend muss das
%          Skript nur noch ausgefuehrt werden.

clear data

%% <Elektroblech Name>
%{%
% Quelle: [xxx]

% Bezeichnung
data.Bezeichnung = '<Elektroblech Name>';

% Dichte des Blechwerkstoffs rho_Fe [kg/m^3]
data.rho_Fe = 0;

% electric steel sheet
% Feldstaerke H [A/m]
data.H = 0;

% Magnetische Induktion [T]
data.B = 0;

% Relative Permeabilitaet [-]
data.mu_r = 0;

% Frequenz (der Eisenverluste) [Hz]
data.p_vFe_f_vec = 0;

% Magnetische Induktion (der Eisenverluste) [T]
data.p_vFe_B_vec = 0;

% Spezifische Eisenverluste [W/kg]
data.p_vFe = 0;

%}

%% VACOFLUX 48
%{
% Quelle: [Datenblatt vacuumschmelze VACOFLUX 48 - Strip 0.35mm]
% V/A: BH-Kurve nur fuer DC vermessen -> keine Frequenzabhaengigkeit

% Bezeichnung
data.Bezeichnung = 'VACOFLUX 48';

% Dichte des Blechwerkstoffs rho_Fe [kg/m^3]
data.rho_Fe = 8120;

% electric steel sheet
% Feldstaerke H [A/m]
data.H = [20 25 30 35.0131233595801 40 45.0131233595801 50 54.9868766404200 60 64.9868766404200 70 80 90 100 120 140 160 180 200 250 300 349.868766404199 400 449.868766404200 500 549.868766404199 599.737532808399 649.606299212598 699.475065616798 799.737532808399 899.212598425197 999.475065616798 1402.09973753281 1601.31233595801 1800.52493438320 2001.04986876640 2199.73753280840 2401.83727034121 2599.21259842520 2797.90026246719 3000 3503.93700787402 4018.37270341207 5498.68766404200 6524.93438320210 7976.37795275591 13976.3779527559 14971.1286089239 16005.2493438320];

% Magnetische Induktion [T]
data.B = [0.0906633897820952 0.157515889990914 0.273117086758649 0.485303231067704 0.853772035038643 1.28962866944155 1.47991707716640 1.57440256518322 1.63682332888747 1.68229943212849 1.71865137414349 1.77310715662577 1.81348577539586 1.84591340355073 1.89421495808777 1.93104623108155 1.95862307237564 1.98333234328390 2.00308853807121 2.04426734620737 2.07071727008909 2.09221448602372 2.11136477165859 2.12712684682158 2.14054199168481 2.15213305537752 2.16176895744353 2.17075313896731 2.17869456762352 2.19236010239691 2.20289958761076 2.21226376680542 2.23915819271357 2.24641995129824 2.25250861290689 2.25742049580090 2.26116443615296 2.26463811576773 2.26799470553307 2.27095658123445 2.27325716387320 2.27744817290733 2.28082766338184 2.28736287911407 2.29021890028695 2.29279444690692 2.30307311250525 2.30471405432047 2.30624454397718];

% Relative Permeabilitaet [-]
data.mu_r = [3609.61419421139 5016.11903540410 7246.88826730563 11032.1551331482 16987.4874901428 22801.2256025367 23555.8443536036 22787.1192173128 21711.2761806039 20602.2585872578 19540.2226515048 17639.6555719274 16036.9672772834 14691.5446686236 12563.6356113177 10978.5022844467 9743.62452209001 8770.48659576905 7972.26859100826 6509.33759689454 5494.98068856337 4760.96365586748 4202.65928488391 3764.91598223267 3409.01092261355 3116.81763930181 2870.62242595304 2661.42851703854 2480.87707012450 2183.72916884556 1951.72901068130 1763.62071665966 1273.08754225773 1118.59450613230 997.769649913936 899.960362044889 820.228980715309 752.550540212076 696.601606489636 648.134695412233 605.232724805272 519.460916619058 453.914127112172 333.261618553199 281.545426047837 230.976439883044 133.362885185159 124.737203958171 116.898107435854];

% Frequenz (der Eisenverluste) [Hz]
data.p_vFe_f_vec = [50 60 100 200 300 400 500 600 700 800 900 1000];

% Magnetische Induktion (der Eisenverluste) [T]
data.p_vFe_B_vec = [0.5 1.0 1.5 2.0];

% Spezifische Eisenverluste [W/kg]
data.p_vFe = [0.262	0.325	0.592	1.419	2.459	3.707	5.104	6.666	8.363	10.44	12.4	14.59; ...
              0.75	0.93	1.768	4.541	8.057	12.25	17.19	22.79	29.15	36.07	43.88	52.36; ...
              1.404	1.769	3.488	9.527	17.98	28.83	41.97	57.9	75.96	96.37	119.8	144.6; ...
              2.25	2.865	5.903	17.62	35.49	59.56	89.67	125.2	167.4	214.8	269.1	334.7];

%}

%% VACOFLUX 50
%{
% Quelle: [Datenblatt vacuumschmelze VACOFLUX 50 - Strip 0.35mm]
% V/A: BH-Kurve nur fuer DC vermessen -> keine Frequenzabhaengigkeit

% Bezeichnung
data.Bezeichnung = 'VACOFLUX 50';

% Dichte des Blechwerkstoffs rho_Fe [kg/m^3]
data.rho_Fe = 8120;

% Feldstaerke H [A/m]
data.H = [9.96325459317585 20 25 30 35.0131233595801 40 45.0131233595801 50 55.0131233595801 60 65.0131233595801 70 80 90 100 120 140 160 180 200 250 300 348.818897637795 398.687664041995 449.081364829396 498.950131233596 548.818897637795 598.950131233596 648.818897637795 698.687664041995 798.687664041995 898.950131233596 998.687664041995 1198.68766404199 1398.68766404199 1598.68766404199 1798.68766404199 1998.68766404199 2198.68766404199 2398.68766404199 2598.68766404199 2797.90026246719 2997.37532808399 3498.68766404199 3997.37532808399 4498.68766404200 4997.37532808399 5997.37532808399 6498.68766404200 6997.37532808399 7997.37532808399 9994.75065616798 10994.7506561680 11992.1259842520 15992.1259842520];

% Magnetische Induktion [T]
data.B = [0.0227906435679169 0.0598563562670374 0.0855878928691798 0.119616159583572 0.169383777954682 0.250386211802010 0.389099657405916 0.642747284030209 1.02047976330908 1.31817326597046 1.44078129429765 1.51543916982385 1.60997203316966 1.67192802975120 1.71741257684519 1.78446756659198 1.83322094579703 1.86952922983371 1.89949295554926 1.92530831620869 1.97605950458248 2.00887511462538 2.03522716243537 2.05864827982520 2.07694360989660 2.09426419043919 2.10816847034732 2.12121800691347 2.13219402913492 2.14170592251304 2.15865485855784 2.17194280431093 2.18327991463964 2.20168242313724 2.21483846994622 2.22482223759463 2.23273182271498 2.23917727899199 2.24440262789956 2.24865189091156 2.25192506802801 2.25495622822068 2.25786470949303 2.26281102683825 2.26690995085856 2.26990409641266 2.27217285937881 2.27633767126655 2.27786768797731 2.27928237578484 2.28186104809229 2.28568295643453 2.28740755358337 2.28889481109204 2.29432907084407];

% Relative Permeabilitaet [-]
data.mu_r = [1822.33646165332 2383.63461011837 2726.37311006614 3174.94304484198 3851.76318679601 4983.30127756926 6880.81245222581 10231.6666075061 14763.4490265561 17484.8417938418 17637.5010636096 17229.8569257436 16016.7138236925 14785.1153345428 13668.7609129313 11835.6439500378 10422.2464128190 9300.30143448934 8399.62845638453 7662.58425378988 6292.01862623669 5330.73260872500 4645.07223518064 4111.05756417316 3682.38167622780 3342.16426071525 3058.82141732605 2820.31000689819 2617.15637960856 2441.33531568889 2152.80774658053 1924.68776602255 1741.70786636119 1463.66033233577 1262.14454646589 1109.47000492274 989.830234114968 893.551187671192 814.346275310537 748.025583322442 691.614328576752 643.377230711055 601.467527755225 516.701074219904 453.309389352709 403.550137650104 363.843338837406 304.066653044956 280.954413848105 261.236670184707 229.081275303628 184.010266573881 167.583153750643 153.912580892586 116.192617170969];

% Frequenz (der Eisenverluste) [Hz]
data.p_vFe_f_vec = [50 60 100 200 300 400 500 600 700 800 900 1000];

% Magnetische Induktion (der Eisenverluste) [T]
data.p_vFe_B_vec = [0.5 1.0 1.5 2.0];

% Spezifische Eisenverluste [W/kg]
data.p_vFe = [0.289439	0.356055	0.646035	1.51366	2.55985	3.76652	5.1218	6.58498	8.17415	9.96743	11.832	13.8546; ...
              0.83864	1.04285     2.00599     5.01346	8.68793	13.247	18.4751	24.5549	31.3967	39.0839	48.4845	56.8055; ...
              1.566     1.946       3.755       9.791	18.43	29.65	43.54	58.09	76.55	97.84	121.1	147.3; ...
              2.436     3.075       6.146       17.61	34.75	57.76	86.42	121.5	162.1	208.5	263.1	321.9];

%}

%% JFESteel35JN210
%{
% Quelle: [Datenblatt JFESteel35JN210 aus Ansys]
% V/A: BH-Kurve nur fuer DC vermessen -> keine Frequenzabhaengigkeit

% Bezeichnung
data.Bezeichnung = 'JFESteel35JN210';

% Dichte des Blechwerkstoffs rho_Fe [kg/m^3]
data.rho_Fe = 0;

% Feldstaerke H [A/m]
data.H = [0 10 20 30 40 50 60 70 80 100 125 150 175 200 250 300 400 500 800 1000 1500 2000 2500 3000 4000 5000 8000 10000 15000 20000 30000 50000 80000 100000 110000 120000 130000 140000 150000 160000 170000 180000 190000 200000 210000 220000 230000 240000];

% Magnetische Induktion [T]
data.B = [0 0.0400000000000000 0.0750000000000000 0.161000000000000 0.330000000000000 0.522000000000000 0.704000000000000 0.842000000000000 0.944000000000000 1.07600000000000 1.17600000000000 1.23700000000000 1.27200000000000 1.29800000000000 1.33600000000000 1.35900000000000 1.39400000000000 1.41800000000000 1.46400000000000 1.48300000000000 1.52100000000000 1.55200000000000 1.57500000000000 1.59600000000000 1.63700000000000 1.67700000000000 1.77500000000000 1.81900000000000 1.88800000000000 1.93100000000000 1.97900000000000 2.03500000000000 2.08100000000000 2.11500000000000 2.12500000000000 2.14500000000000 2.16500000000000 2.18000000000000 2.19000000000000 2.20500000000000 2.21500000000000 2.23000000000000 2.25000000000000 2.26500000000000 2.28500000000000 2.29800000000000 2.31000000000000 2.33000000000000];  

% Relative Permeabilitaet [-]
data.mu_r = 0;

% Frequenz (der Eisenverluste) [Hz]
data.p_vFe_f_vec = 0;

% Magnetische Induktion (der Eisenverluste) [T]
data.p_vFe_B_vec = 0;

% Spezifische Eisenverluste [W/kg]
data.p_vFe = 0;

%}

%% M800-50A
%{
% Quelle: [Pyr14, S.324], [https://cogent-power.com/cms-data/downloads/m800%2D50a%2Epdf]
% V/A: BH-Kurve nur fuer DC vermessen -> keine Frequenzabhaengigkeit

% Magnetische Feldkonstante mu_0 [H/m]
const.mu_0 = pi * 4e-7;

% Bezeichnung
data.Bezeichnung = 'M800-50A';

% Dichte des Blechwerkstoffs rho_Fe [kg/m^3]
data.rho_Fe = 7800;

% Magnetische Induktion [T]
data.B = 0:0.1:2.0;

% Feldstaerke H [A/m]
data.H = [1835.2.*data.B(data.B<1.5).^5 - 6232.3.*data.B(data.B<1.5).^4 + 7806.7.*data.B(data.B<1.5).^3 - 4376.3.*data.B(data.B<1.5).^2 + 1227.6.*data.B(data.B<1.5) 0.011637.*exp(7.362.*data.B(data.B>=1.5))];

% Relative Permeabilitaet [-]
data.mu_r = (data.B ./ data.H) ./ const.mu_0;

% Frequenz (der Eisenverluste) [Hz]
data.p_vFe_f_vec = [50 100 200 400 1000 2500];

% Magnetische Induktion (der Eisenverluste) [T]
data.p_vFe_B_vec = 0.1:0.1:1.8;

% Spezifische Eisenverluste [W/kg]
data.p_vFe = [0.05 0.18 0.43 0.7 1.01 1.35 1.72 2.13 2.56 3.05 3.59 4.20 4.91 5.70 6.60 7.54 8.30 8.83;
              0.05 0.18 0.43 0.7 1.01 1.35 1.72 2.13 2.56 3.05 3.59 4.20 4.91 5.70 6.60 7.54 8.30 8.83;
              0.05 0.18 0.43 0.7 1.01 1.35 1.72 2.13 2.56 3.05 3.59 4.20 4.91 5.70 6.60 7.54 8.30 8.83;
              0.05 0.18 0.43 0.7 1.01 1.35 1.72 2.13 2.56 3.05 3.59 4.20 4.91 5.70 6.60 7.54 8.30 8.83;
              0.05 0.18 0.43 0.7 1.01 1.35 1.72 2.13 2.56 3.05 3.59 4.20 4.91 5.70 6.60 7.54 8.30 8.83;
              0.05 0.18 0.43 0.7 1.01 1.35 1.72 2.13 2.56 3.05 3.59 4.20 4.91 5.70 6.60 7.54 8.30 8.83]';

%}

%% M250-35A
%{
% Quelle: [https://cogent-power.com/cms-data/downloads/m250%2D35a%5F1%2Epdf]
% V/A: BH-Kurve nur fuer DC vermessen -> keine Frequenzabhaengigkeit

% Magnetische Feldkonstante mu_0 [H/m]
const.mu_0 = pi * 4e-7;

% Bezeichnung
data.Bezeichnung = 'M250-35A';

% Dichte des Blechwerkstoffs rho_Fe [kg/m^3]
data.rho_Fe = 7600;

% Magnetische Induktion [T]
data.B = 0.1:0.1:1.8;

% Feldstaerke H [A/m]
data.H = [26.8 35.7 41.8 47.5 53.4 60.0 67.9 77.5 90.0 107.0 133.0 179.0 284.0 642.0 1810.0 4030.0 7290.0 11700.0];

% Relative Permeabilitaet [-]
data.mu_r = (data.B ./ data.H) ./ const.mu_0;

% Frequenz (der Eisenverluste) [Hz]
data.p_vFe_f_vec = [50 100 200 400 1000 2500];

% Magnetische Induktion (der Eisenverluste) [T]
data.p_vFe_B_vec = 0.1:0.1:1.8;

% Spezifische Eisenverluste [W/kg]
data.p_vFe = [0.02 0.06 0.13 0.21 0.31 0.41 0.52 0.66 0.81 0.98 1.15 1.37 1.65 2.00 2.35 2.65 2.87 3.06; ...
              0.04 0.14 0.31 0.51 0.75 1.01 1.31 1.64 2.00 2.41 2.87 3.40 4.03 4.83 5.72 NaN  NaN  NaN; ...
              0.08 0.33 0.73 1.23 1.82 2.49 3.26 4.12 5.07 6.14 7.33 8.69 10.3 12.4 14.7 NaN  NaN  NaN; ...
              0.21 0.90 1.93 3.24 4.81 6.69 8.82 11.2 14.0 17.1 20.6 24.6 29.2 35.1 41.6 NaN  NaN  NaN; ...
              0.98 3.65 7.58 12.7 18.8 26.3 35.2 45.7 58.1 72.6 89.6 NaN  NaN  NaN  NaN  NaN  NaN  NaN; ...
              4.09 14.8 30.6 51.7 78.8 113  155  208  273  352  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN]';

%}

%% Datei erzeugen
path = ['5_Materialien/1_Elektroblech/', data.Bezeichnung];
save(path,'data');
