function [RES, Optimierer]=Powertrain_Calculation(folder_id, file_id, LDS_values, Maschinendaten, vehicle, CycEff_export)

%load  (['3_Ergebnisse/',folder_id,'/','1_Entwurf/', 'Entwurf_', file_id,'.mat']); %,'motor');

Achsmoment_Motor=               Maschinendaten.Entwurf.Bemessungswerte.M_N;     % Motor torque installed on the axis % Nominal machine torque in Nm
Nenndrehzahl_Motor=             Maschinendaten.Entwurf.Bemessungswerte.n_N;     % Rated speed of the motors in 1/min
Maximale_Drehzahl=              Maschinendaten.Entwurf.Bemessungswerte.nmax;     % Maximum speed of the e-machine(s) in 1/min
Maschinentyp=                   vehicle.Maschinendaten.type;     % PSM oder ASM
n_Gaenge=                       1;     % number of gears
i_Gaenge=                       vehicle.GEARBOX{1, 1}.gear_ratio;     % gear ratio
az=                             1;     % 1 for axle center drive
% etv=                            varargin{8};     % 1 for electric torque vectoring (only in conjunction with axle center drive and without differential)
% TS=                             varargin{9};     % 1 for Torque Splitter (only in combination with axle center drive and without differential)
% OD=                             varargin{10};    % 1 for open differential (only in conjunction with axle-center drive without TV)
Stkzahl=                        100000;    % Number of units produced per year (maximum 2 million) %100000 as assumption
Standort=                       1;    % Germany=1, USA=2, Romania=3, Czech Republic=4
Magnetpreise=                   3;    % stable=1, decreasing=2, increasing=3 pcs
%Optimierung                 =    varargin{14};    % Optimization of the overall topology


% This script serves as a tool for calculating the manufacturing costs, mass, efficiency and
% geometry of electric powertrains. The results
% are calculated for a selected powertrain topology of an axle.
% The powertrain topology
% can be freely modified in the "Data_Initialization" script. 
% The results of the efficiency map calculation are available under
% "Results.Map" and correspond to the 4 outputs of the
% "Pesce model". However, the maps are calculated with the models of Horlbeck
% models.

% All results are stored in the structure-array "RES".


RES.em.typ_EM =             Maschinentyp;
RES.trans.n_gears=          n_Gaenge;
RES.trans.i_gears =         i_Gaenge;
RES.em.nmax=                Maximale_Drehzahl; % in 1/min
RES.em.nnenn=               Nenndrehzahl_Motor; % in 1/min
RES.em.Mnenn_achs=          Achsmoment_Motor; % in Nm
RES.em.U_RMS=               Maschinendaten.Entwurf.Bemessungswerte.U_N;                            %V 
RES.em.U_linetoline=        Maschinendaten.Entwurf.Bemessungswerte.U_N;                    %V linetoline voltage
RES.akt.az=                 az;
% RES.akt.etv=                etv;
% RES.akt.TS=                 TS;
% RES.akt.OD=                 OD;
RES.cost.Stkzahl=           Stkzahl;
RES.cost.Standort=          Standort;
RES.cost.Magnetpreise=      Magnetpreise;
%RES.trans.etv.delta_M=      delta_M_TV_max;

% Halve M_nenn_mot for drive close to wheel
if (RES.akt.az == 0)
    RES.em.Mnenn_mot=Achsmoment_Motor*0.5; % in Nm
else
    RES.em.Mnenn_mot=Achsmoment_Motor; % in Nm
end
%% Motor
disp('Auslegung EM')
tic

%% Machine geometry
D = Maschinendaten.Entwurf.Geometrie.D_2a; % Outside diameter of armature in m
D_rot = D; % Rotor diameter in m
D_a = Maschinendaten.Entwurf.Geometrie.D_1a; % Stator outer diameter in m
A_nuten_Stator = Maschinendaten.Entwurf.Geometrie.Nut_1.A_1n; % Stator groove cross section in mm^2
N_S = Maschinendaten.Entwurf.Wicklung.N_1; % Stator number of slots, [-]
D_i = Maschinendaten.Entwurf.Geometrie.D_2i;% Rotor inner diameter in m
P_n = Maschinendaten.Entwurf.Bemessungswerte.P_N/1000; % rated power
P_n_kW = P_n; % rated power in kW
M_n = Maschinendaten.Entwurf.Bemessungswerte.M_N; %Rated Torque in Nm

A_1S = Maschinendaten.Entwurf.Wicklung.A_1L; % Conductor cross sectional area in mm^2
z_nS = Maschinendaten.Entwurf.Wicklung.z_1n; % Number of conductors per slot stator, [-]
S_S_max = Maschinendaten.Entwurf.Wicklung.S_1; % maximum current density in the stator in A/mm^2
if strcmp(Maschinentyp,'ASM')
    cphi=Maschinendaten.Entwurf.EMAG.cos_phi_N; %[-]
else
    cphi=Maschinendaten.Entwurf.Bemessungswerte.cos_phi_N; %[-]
end

A_wS = A_1S*N_S*z_nS; % Winding cross section stator in mm^2
          
if strcmp(Maschinentyp,'ASM')
A_nuten_Rotor = Maschinendaten.Entwurf.Geometrie.Nut_2.A_2n; % Rotor slot cross section, [mm^2]
N_L = Maschinendaten.Entwurf.Wicklung.N_2; % Rotor number of slots, [-]
else
end

if strcmp(Maschinentyp,'ASM')
S_L_max = Maschinendaten.Entwurf.Richtwerte.S_2; % maximum current density S_2 in the rotor in A/mm^2
A_wL = A_wS*(S_S_max/S_L_max)*cphi; % Winding cross section rotor in mm^2
A_nL = Maschinendaten.Entwurf.Geometrie.Nut_2.A_2n; %in mm^2
A_nuten_Rotor=A_nL*N_L; %/1000000; in mm^2
else
end

t_gehaeuse=0.01; % in m
L_gehaeuseueberhang=0.0459; % in m
l_WK=0.03; % length of the winding overhead in m
l_fe = Maschinendaten.Entwurf.Geometrie.l_Fe; % Sheet bundle length in m
L_M = l_fe+2*l_WK+L_gehaeuseueberhang; % Length with housing in m
MV = (D_a^2*pi/4)*L_M; % aktive Volume in m^3
t_luft=0.005; % in m, Air between machine and housing
D_aussen_EM=D_a+2*(t_gehaeuse+t_luft); % Housing outer diameter in m
D_M = D_aussen_EM; %in m
V_Gehaeuse=((D_aussen_EM^2-(D_aussen_EM-2*t_gehaeuse)^2)*pi/4)*(l_fe+2*l_WK)+D_aussen_EM^2*pi/4*2*t_gehaeuse; % Volume of the housing in m^3, rear part represents lid
A_nS = Maschinendaten.Entwurf.Geometrie.Nut_1.A_1n; % in mm^2

A_Luftspalt = ((D+(2*Maschinendaten.Entwurf.Geometrie.delta_i/1000))^2*pi/4)-(D^2*pi/4); % AIR GAP ONLY! in m^2

if strcmp(Maschinentyp,'PMSM')
A_Magnete = ((Maschinendaten.Entwurf.Geometrie.b_PM * Maschinendaten.Entwurf.Geometrie.h_PM)/1000000)*Maschinendaten.Entwurf.Geometrie.anzahl_PM; % Width x height in m^2
end

[d_welle_a, d_welle_i]=Rotor_Design(M_n, D_i); % in m

A_Welle=(d_welle_a^2-d_welle_i^2)*pi/4; % Cross-sectional area of the rotor shaft in m^2

typ_EM = Maschinendaten.Entwurf.Optionen.Maschinentyp;
                  

% [L_M, D, MV, D_a, l_fe, A_wS, A_wL, l_WK, D_rot, D_M, D_i, P_n_kW, V_Gehaeuse, A_nuten_Stator, A_nuten_Rotor, A_Luftspalt_und_Magnete, A_Welle]    =...
%     Maschinegeometry(RES.em.Mnenn_mot, RES.em.nnenn, RES.em.typ_EM, RES.em.U_linetoline, Maschinendaten);

RES.em.Motoraussendurchmesser_in_Meter=D_M; %in m
RES.em.Maschinenlaenge_in_Meter=L_M; %in m
RES.em.Statoraussendurchmesser_in_Meter=D_a; %in m
RES.em.Rotoraussendurchmesser_in_Meter=D_rot; %in m

% Mass of the individual components, all in kg
if strcmp(Maschinentyp,'ASM')
A_Magnete=0;
    [RES.em.Motormasse, RES.em.Masse_Aluminium, RES.em.Masse_Elektroblech, ...
    RES.em.Masse_Kupferdraht, RES.em.Masse_Magnete, RES.em.Masse_Rotorkaefig,...
    RES.em.Masse_Stahl] = MassCalc_EM(D, D_i, MV, A_wS, A_wL, l_WK, D_a, l_fe, D_rot, A_nuten_Stator, A_nuten_Rotor, A_Luftspalt, A_Magnete, A_Welle, RES.em.typ_EM, V_Gehaeuse, Maschinendaten);
else
    A_wL = 0;
    A_nuten_Rotor = 0;
    [RES.em.Motormasse, RES.em.Masse_Aluminium, RES.em.Masse_Elektroblech, ...
    RES.em.Masse_Kupferdraht, RES.em.Masse_Magnete, RES.em.Masse_Rotorkaefig,...
    RES.em.Masse_Stahl] = MassCalc_EM(D, D_i, MV, A_wS, A_wL, l_WK, D_a, l_fe, D_rot, A_nuten_Stator, A_nuten_Rotor, A_Luftspalt, A_Magnete, A_Welle, RES.em.typ_EM, V_Gehaeuse, Maschinendaten);
end

% Costs in Euro
if strcmp(Maschinentyp,'ASM')
    [RES.em.K_ges, MK, K_Anbauteile, KO_FT] = Cost_Calculation_Motor(MV, A_wS, l_WK, typ_EM, D_a, l_fe, D, A_wL, P_n, Stkzahl, Magnetpreise, az, Standort, V_Gehaeuse, A_nuten_Stator, A_nuten_Rotor, A_Luftspalt, A_Magnete, D_i, Maschinendaten);            %Berechnung der Herstellkosten des Motors an der VA
else
    A_wL = 0;
    A_nuten_Rotor = 0;
    [RES.em.K_ges, MK, K_Anbauteile, KO_FT] = Cost_Calculation_Motor(MV, A_wS, l_WK, typ_EM, D_a, l_fe, D, A_wL, P_n, Stkzahl, Magnetpreise, az, Standort, V_Gehaeuse, A_nuten_Stator, A_nuten_Rotor, A_Luftspalt, A_Magnete, D_i, Maschinendaten);            %Berechnung der Herstellkosten des Motors an der VA
end

RES.cost.Materialkosten = MK;
RES.cost.Anbauteile = K_Anbauteile;
RES.cost.Fertigung = KO_FT;

[RES.em.Jx, RES.em.Jy, RES.em.Jz] = Inertness_Axle_EM(RES.em.typ_EM, RES.em.Mnenn_mot, RES.em.nnenn);
% [RES.trans.Jx_em_Steuer, RES.trans.Jy_em_Steuer, RES.trans.Jz_em_Steuer] =  Inertness_Axle_EM(RES.trans.typ_em_Steuer, RES.trans.M_nenn_em_Steuer,  RES.trans.n_nenn_em_Steuer);

RES.em.Jred = inertness_EM(RES.em.Mnenn_mot,RES.em.nnenn,RES.em.typ_EM,RES.trans.i_gears);


RES.em.M_max_mot = LDS_values.T_max(1,1);                           % Maximum moment, first element from full load curve vector
LDS_values.eff_n_axis(1,1) = 0;
LDS_values.eff_n_axis = LDS_values.eff_n_axis';
RES.em.MaxMotorTrqCurve_w = LDS_values.eff_n_axis;                  % Axis speed motor range incl.0, first line from 0 high
RES.em.MaxMotorTrqCurve_M = LDS_values.T_max;                       % Full load curve motor range, first line
RES.em.MaxGeneratorTrqCurve_w = LDS_values.eff_n_axis;              % Axis speed generator range,like engine range from 0 up, first line
RES.em.MaxGeneratorTrqCurve_M = LDS_values.T_max;                   % Full load curve generator range,like engine range first line
RES.em.MotorEffMap3D_w = LDS_values.eff_n_axis;                     % Axis speed generator range,like engine range from 0 up, first line
ticker = Maschinendaten.Analyse.Optionen.M_tics;
RES.em.MotorEffMap3D_M = LDS_values.eff_T_axis(1,ticker:end);       % Axis torque from 0 to high, first line
LDS_values.eff_recu(isnan(LDS_values.eff_recu))=1.0e-3;
RES.em.MotorEffMap3D = LDS_values.eff_recu;                         % Map only motor area folded down, NaN replaced with 1e^-3
RES.em.GeneratorEffMap3D_w = LDS_values.eff_n_axis;                 % Axis speed generator range,like engine range from 0 up, first line 
RES.em.GeneratorEffMap3D_M = LDS_values.eff_T_axis(1,ticker:end);   % Axis torque from 0 to high, first line 
RES.em.GeneratorEffMap3D = LDS_values.eff_recu;                     % Map only motor area folded down, NaN replaced with 1e^-3 

% Calculate maximum axis torque
if az==1
    RES.em.M_max_achse=RES.em.M_max_mot;
else
    RES.em.M_max_achse=RES.em.M_max_mot*2;
end
toc

%% Total costs
%Kosten_gesamt=RES.em.K_ges + RES.trans.K_ges + RES.LE.K_ges;
if az == 1
RES.cost.Gesamtkosten = RES.em.K_ges; % + RES.trans.K_ges + RES.inv.Kosten; % Costs in Euro
else
RES.cost.Gesamtkosten = (RES.em.K_ges)*2; % + RES.trans.K_ges + RES.inv.Kosten) * 2; % Costs in Euro
end

% Data
RES.Masse.Topologie_Gesamtmasse = RES.em.Motormasse; %+RES.trans.masse+RES.trans.etv.m_em_Steuer; % in kg


%% Centers of gravity [m] relative to the respective axis
%coordinates (axis-specific KOSY, [m]):

%length of motors from regressions
[l_EM] = RES.em.Maschinenlaenge_in_Meter; %m
 disp(['Cost of machine: ' char(9) char(9) char(9) char(9) num2str(round(RES.em.K_ges,2)) char(9) 'Euro']);

folder_id = [Maschinendaten.Entwurf.Optionen.folder_id];
file_id = [Maschinendaten.Entwurf.Optionen.file_id];

% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.cost.Stkzahl,'Kosten','C3:C3')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.cost.Standort,'Kosten','C4:C4')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.cost.Magnetpreise,'Kosten','C5:C5')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.cost.Gesamtkosten,'Kosten','C6:C6')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.em.Motoraussendurchmesser_in_Meter,'Kosten','C7:C7')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.em.Maschinenlaenge_in_Meter,'Kosten','C8:C8')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.em.Motormasse,'Kosten','C9:C9')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.em.Masse_Aluminium,'Kosten','C10:C10')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.em.Masse_Elektroblech,'Kosten','C11:C11')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.em.Masse_Kupferdraht,'Kosten','C12:C12')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.em.Masse_Magnete,'Kosten','C13:C13')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.em.Masse_Rotorkaefig,'Kosten','C14:C14')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.em.Masse_Stahl,'Kosten','C15:C15')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.em.K_ges,'Kosten','C16:C16')

%% Parameter for optimization
r_Zylinder = 1/2 * D_M;
V_Zylinder = 2*L_M*(r_Zylinder)^2; %in m^3
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], RES.em.K_ges,'Optimierung','C3:C3')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], V_Zylinder,'Optimierung','C4:C4') % in m
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], (table2array(CycEff_export)*100)','Optimierung','C5:C5')
Zykluseffizienz = table2array(CycEff_export)*100;

Optimierer.Kosten_in_Euro = RES.em.K_ges;
Optimierer.Volumen_in_m3 = V_Zylinder;
if isnan(Zykluseffizienz)
    Optimierer.Zykluseffizienz_in_Prozent = -Inf;
else
    Optimierer.Zykluseffizienz_in_Prozent = Zykluseffizienz;
end
%Optimierer.Zykluseffizienz = table2array(CycEff_export)*100;

disp(['All done!']); %neu
%[l_GTR] = RES.trans.Getriebelaenge_in_Meter; %m
%[b_GTR] = RES.trans.Getriebebreite_in_Meter; %m

%Schwerpunktberechnung EM und Getriebe
% [RES.em.CoG_x, RES.em.CoG_y, RES.em.CoG_z] = SP_EM(l_EM,l_GTR,b_GTR,RES.akt.az); %m
% [RES.trans.CoG_x, RES.trans.CoG_y, RES.trans.CoG_z] = SP_GTR(l_EM, l_GTR, b_GTR, RES.akt.az);
%[RES.trans.CoG_x_AW, RES.trans.CoG_y_AW, RES.trans.CoG_z_AW] = SP_AW(b_GTR,RES.akt.az,config.segment_parameter.Breite,b_Reifen);

