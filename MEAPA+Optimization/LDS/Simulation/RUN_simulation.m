function [FZG_LDS1, Optimierer] = RUN_simulation(LDS_values, handles_GUI_LDS,FZD_LDS, Maschinendaten) %FZG_LDS1

% This script runs the simulation of the vehicle maximum acceleration for
% the practical course "Praktikum Fahrzeugkonzeptentwicklung". Follow the
% steps to determine acceleration time 0-100km/h and maximum velocity.

%% Matlab Tip:
% Use "Crtl + d" or right click and "Open function" to access a function 
% within a script - without searching it in the folder structure

%% 1) Initialization -------------------------------------------------------------------------------
% Specifiy the powertrain architecture. Use a string to create a unique
% name for a powertrain archticture. Specify the number of gearboxes and
% motors at each axle. For example:
%   - GM_X:     one gearbox, one motor at the front axle
%   - X_GM:     one gearbox, one motor at the rear axle 
%   - GM_GM:    one gearbox and one motor at each axle
%   - 2G2M_X    two gearboxes, two motors at the front axle
%IMPORTANT NOTICE: For the simulation there is no difference between "front" and "rear" axle. GM_X
% and X_GM are the same for the simulation. Only Important for Package  For this longitudinal dynamic simulation, both architectures
%       deliver the same results, because no weight distribution and no
%       wheel slip are considered.



architecture = 'GM_X'; %handles_GUI_LDS.AnzMasch.String; %'GM_X'; %Choose the powertrain architecture

%Choose which simulation you want to run: 
whichSim = 'E'; % 'A' for acceleration, 'E' for energy consumption

%% Initialize vehicle
[vehicle] = Para_vehicle(architecture, LDS_values, FZD_LDS);   %Open this function and define all vehicle parameters for your concept
%load('C:\00_All\10_Praktikum_FKE\Vehicle\Saved_Vehicles\vehicle_eGolf.mat')

%% Initialize driving cycle (for energy consumption simulation)
name_cycle = 'WLTP_class_3.mat';
load(name_cycle);               %load different driving cycle if you want
driving_cycle = dc;             %rename duty cycle

%% Initialize simulation parameters for Acceleration Simulation
v_max_sim = 350;    %maximum simulated velocity in km/h. When vehicle reaches this velocity, the simulation ends.
t_sim = 100;        %simulated time in s. Simulation ends after this time.


%% Initialize simulation parameters for Consumption Simulation
if (exist('handles_GUI_LDS.visual_LDS'))
    visualize_ECS = handles_GUI_LDS.visual_LDS;    %binary to determine if the results of energy consumption simulation should be visualized (1 = yes, 0 = no)
else
    visualize_ECS = 1;    %binary to determine if the results of energy consumption simulation should be visualized (1 = yes, 0 = no)
end

%% 2) Simulation -----------------------------------------------------------------------------------
if whichSim == 'A'
%% 2.1) Acceleration Simulation 
    [RESacc, v_acc] = Simulation_acceleration(vehicle, v_max_sim, t_sim); %function that simulates the maximum acceleration

elseif whichSim == 'E'
%% 2.1) Energy Consumption Simulation 
tic
    [P_el_min, E_EL_MIN, data_R, F_Ges, v, V] = Simulation_energy_consumption(driving_cycle, vehicle, architecture); %function that simulates the energy consumption of a driving cycle

toc
else
    disp(['Please select a valid simulation type in this script. Choose "A" for acceleration' ...
        ' and "E" for energy consumption.']);
end

%speed-longitudinal force plot
%figure; plotyy(dc.time,v,dc.time,F_Ges), title('speed-longitudinal force'); xlabel('time in s'); ylabel('vehicle speed in m/s'); ylabel('total driving resistance in N');
%figure; x = dc.time; y = v; yyaxis left; plot(x,y); title('speed-longitudinal force'); xlabel('time in s'); ylabel('vehicle speed in m/s'); x2=dc.time; y2=F_Ges; yyaxis right; plot(x2,y2);title('speed-longitudinal force'); xlabel('time in s'); ylabel('total driving resistance in N')
%% 3) Result and Visualization ---------------------------------------------------------------------
if whichSim == 'A'
%% 3.1) Acceleration Simulation 
    disp([newline,newline, 'Determine the 0-100km/h time and the maximum' ...
        'speed from the figure or calculate it using the result vector "v_acc".']);
    
    plot(v_acc.t,v_acc.v)
    xlabel('time in s')
    ylabel('velocity in km/h')
    title('Acceleration simulation')
    
%     plot(RESacc.v_acc, RESacc.F_drive, RESacc.F_aero, RESacc.F_fric, RESacc.F_acc)
%     xlabel('time in s')
%     ylabel('velocity in km/h')
%     title('Acceleration simulation')
    
FZG_LDS1 = 'Please Wait...';
disp(FZG_LDS1);
    
elseif whichSim == 'E'
%% 3.2) Energy Consumption Simulation 
    Results.P = cumsum(P_el_min); % Required Power in Ws
    Results.E = cumsum(E_EL_MIN);  % Required Energy in W
    Results.S = cumsum(driving_cycle.speed); % Distance in m
    Results.EC = Results.E(end)/(1000*3600)*100/(Results.S(end)/1000); % Energy consumption in kWh/100km; 
    Results.Range = vehicle.battery_cap/Results.EC*100; %Range of vehicle in km
    Cycle = string(name_cycle); %Zyklusname
 
    [meanEffMix, effLoadPointsMix] = meanEfficiency(vehicle, data_R.Tn); %calculate the mean efficiency of the load points %neu
    
    disp([newline, 'Energy consumption simulation was successful', newline ...
        'Energy consumption cycle: ' char(9) num2str(round(Results.E(end)/1000/3600,2)) char(9) ' kWh', newline ...
        'Distance driven in cycle:' char(9) num2str(round(Results.S(end),2)) ' m', newline ...
        'Consumption : ' char(9) char(9) char(9) char(9) num2str(Results.EC) char(9)  ' kWh/100km' , newline ...
        'Range in cycle: ' char(9) char(9) char(9) num2str(round(Results.Range,2)) char(9) ' km', newline...
        'Cycle Efficiency:', char(9) char(9) char(9) num2str(meanEffMix*100) char(9) ' %']); %neu
    
    disp(['Cycle name: ' char(9) char(9) char(9) char(9) num2str(Cycle) '']); %neu

    
%% Sum of Nominal and overload load points

%noch ersetzen:
Maschinendaten.Entwurf.Bemessungswerte.nmax = Maschinendaten.Analyse.Optionen.n_max; % Konfiguration für Auslegung mit GUI
Mn_lp = Maschinendaten.Entwurf.Bemessungswerte.M_N; %Nm
neck_lp = Maschinendaten.Entwurf.Bemessungswerte.n_N; % nominal velocity load points
nmax_lp = Maschinendaten.Entwurf.Bemessungswerte.nmax;

CountDatapoints_all = 0;
CountDatapoints_overload1 = 0;
CountDatapoints_overload2 = 0;
CountDatapoints_nominal1 = 0;
CountDatapoints_nominal2 = 0;

CategorizeDataTn = data_R.Tn;
CategorizeDataTn(2,:) = CategorizeDataTn(2,:)*60;
CategorizeData = struct;
CategorizeData.n = data_R.Tn(2,:)*60;
CategorizeData.M = data_R.Tn(1,:);

% Find Torque at max. velocity
xIndex = find(LDS_values.eff_n_axis == max(LDS_values.eff_n_axis), 1, 'first'); %160
yIndex = find(~isnan((LDS_values.eff(:,xIndex))), 1, 'first');
maxMValue = abs(LDS_values.eff_T_axis(yIndex));

% Calculate overload line
m = ((neck_lp - nmax_lp)/(Mn_lp - (2/3)*maxMValue))^-1;
b = Mn_lp - m*neck_lp;

% CategorizeDataTn(1,:)= m*CategorizeDataTn(2,:)+b;

% Calculate all points
for i = 1:length(CategorizeDataTn) 
    if ((CategorizeDataTn(1,i) >= 0) || (CategorizeDataTn(1,i) < 0)) && (CategorizeDataTn(2,i) <= 200000)
        CountDatapoints_all= CountDatapoints_all+1;
    end
end

% Calculate overload points
% Grundstellbereich
for i = 1:length(CategorizeDataTn) 
    if ((CategorizeDataTn(1,i) > Mn_lp) || (CategorizeDataTn(1,i) < Mn_lp*(-1))) && (CategorizeDataTn(2,i) < neck_lp)
        CountDatapoints_overload1= CountDatapoints_overload1+1;
    end
end
% Feldschwächebereich
for i = 1:length(CategorizeDataTn) 
    if ((CategorizeDataTn(1,i) > (m*CategorizeDataTn(2,i) + b))  || (CategorizeDataTn(1,i) < (-m*CategorizeDataTn(2,i) - b))) && (CategorizeDataTn(2,i) >= neck_lp)
        CountDatapoints_overload2= CountDatapoints_overload2+1;
    end
end
CountDatapoints_overload = CountDatapoints_overload1 + CountDatapoints_overload2;

% Calculate nominal points
% Grundstellbereich
for i = 1:length(CategorizeDataTn) 
    if (((CategorizeDataTn(1,i) <= Mn_lp) && (CategorizeDataTn(1,i)>=0)) || ((CategorizeDataTn(1,i) >= Mn_lp*(-1)) && (CategorizeDataTn(1,i) <= 0))) && (CategorizeDataTn(2,i) <= neck_lp)
        CountDatapoints_nominal1= CountDatapoints_nominal1+1;
    end
end
% Feldschwächebereich
for i = 1:length(CategorizeDataTn) 
    if (((CategorizeDataTn(1,i) <= (m*CategorizeDataTn(2,i) + b) && (CategorizeDataTn(1,i) >=0))  || ((CategorizeDataTn(1,i) >= (-m*CategorizeDataTn(2,i) - b)) && (CategorizeDataTn(1,i) <= 0))) && (CategorizeDataTn(2,i) > neck_lp))
        CountDatapoints_nominal2= CountDatapoints_nominal2+1;
    end
end
CountDatapoints_nominal = CountDatapoints_nominal1 + CountDatapoints_nominal2;

Ratio = (CountDatapoints_overload/CountDatapoints_nominal) * 100;

disp(['Overload/Nominal: ' char(9) char(9) char(9) num2str(CountDatapoints_overload) '/' num2str(CountDatapoints_nominal) char(9) '=' char(9) num2str(Ratio) char(9) '%']); %neu

%% Classification operating points for velocity

CountDatapoints_classification1000 = 0;
CountDatapoints_classification2000 = 0;
CountDatapoints_classification3000 = 0;
CountDatapoints_classification4000 = 0;
CountDatapoints_classification5000 = 0;
CountDatapoints_classification6000 = 0;
CountDatapoints_classification7000 = 0;
CountDatapoints_classification8000 = 0;
CountDatapoints_classification9000 = 0;
CountDatapoints_classification10000 = 0;
CountDatapoints_classification11000 = 0;
CountDatapoints_classification12000 = 0;
CountDatapoints_classificationhigher = 0;

for i = 1:length(CategorizeDataTn) 
    if ((CategorizeDataTn(2,i) <= 1000))
        CountDatapoints_classification1000= CountDatapoints_classification1000+1;
    end
    if ((CategorizeDataTn(2,i) > 1000) && (CategorizeDataTn(2,i) <= 2000))
        CountDatapoints_classification2000= CountDatapoints_classification2000+1;
    end
    if ((CategorizeDataTn(2,i) > 2000) && (CategorizeDataTn(2,i) <= 3000))
        CountDatapoints_classification3000= CountDatapoints_classification3000+1;
    end
    if ((CategorizeDataTn(2,i) > 3000) && (CategorizeDataTn(2,i) <= 4000))
        CountDatapoints_classification4000= CountDatapoints_classification4000+1;
    end
    if ((CategorizeDataTn(2,i) > 4000) && (CategorizeDataTn(2,i) <= 5000))
        CountDatapoints_classification5000= CountDatapoints_classification5000+1;
    end
    if ((CategorizeDataTn(2,i) > 5000) && (CategorizeDataTn(2,i) <= 6000))
        CountDatapoints_classification6000= CountDatapoints_classification6000+1;
    end
    if ((CategorizeDataTn(2,i) > 6000) && (CategorizeDataTn(2,i) <= 7000))
        CountDatapoints_classification7000= CountDatapoints_classification7000+1;
    end
    if ((CategorizeDataTn(2,i) > 7000) && (CategorizeDataTn(2,i) <= 8000))
        CountDatapoints_classification8000= CountDatapoints_classification8000+1;
    end
    if ((CategorizeDataTn(2,i) > 8000) && (CategorizeDataTn(2,i) <= 9000))
        CountDatapoints_classification9000= CountDatapoints_classification9000+1;
    end
    if ((CategorizeDataTn(2,i) > 9000) && (CategorizeDataTn(2,i) <= 10000))
        CountDatapoints_classification10000= CountDatapoints_classification10000+1;
    end
    if ((CategorizeDataTn(2,i) > 10000) && (CategorizeDataTn(2,i) <= 11000))
        CountDatapoints_classification11000= CountDatapoints_classification11000+1;
    end
    if ((CategorizeDataTn(2,i) > 11000) && (CategorizeDataTn(2,i) <= 12000))
        CountDatapoints_classification12000= CountDatapoints_classification12000+1;
    end
    if ((CategorizeDataTn(2,i) > 12000))
        CountDatapoints_classificationhigher = CountDatapoints_classificationhigher+1;
    end
end
folder_id = [Maschinendaten.Entwurf.Optionen.folder_id];
file_id = [Maschinendaten.Entwurf.Optionen.file_id];

% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classification1000,'ClassBetr','C3:C3')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classification2000,'ClassBetr','C4:C4')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classification3000,'ClassBetr','C5:C5')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classification4000,'ClassBetr','C6:C6')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classification5000,'ClassBetr','C7:C7')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classification6000,'ClassBetr','C8:C8')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classification7000,'ClassBetr','C9:C9')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classification8000,'ClassBetr','C10:C10')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classification9000,'ClassBetr','C11:C11')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classification10000,'ClassBetr','C12:C12')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classification11000,'ClassBetr','C13:C13')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classification12000,'ClassBetr','C14:C14')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classificationhigher,'ClassBetr','C15:C15')

%% Classification operating points for torque

CountDatapoints_classificationt50 = 0;
CountDatapoints_classificationt100 = 0;
CountDatapoints_classificationt150 = 0;
CountDatapoints_classificationt200 = 0;
CountDatapoints_classificationt250 = 0;
CountDatapoints_classificationt300 = 0;
CountDatapoints_classificationthigher = 0;


for i = 1:length(CategorizeDataTn) 
    if ((CategorizeDataTn(1,i) <= 50) && (CategorizeDataTn(1,i) >= -50))
        CountDatapoints_classificationt50 = CountDatapoints_classificationt50+1;
    end
    if (((CategorizeDataTn(1,i) <= 100) && (CategorizeDataTn(1,i) > 50)) || ((CategorizeDataTn(1,i) >= -100) && (CategorizeDataTn(1,i) < -50))) 
        CountDatapoints_classificationt100 = CountDatapoints_classificationt100+1;
    end
    if (((CategorizeDataTn(1,i) <= 150) && (CategorizeDataTn(1,i) > 100)) || ((CategorizeDataTn(1,i) >= -150) && (CategorizeDataTn(1,i) < -100))) 
        CountDatapoints_classificationt150 = CountDatapoints_classificationt150+1;
    end
    if (((CategorizeDataTn(1,i) <= 200) && (CategorizeDataTn(1,i) > 150)) || ((CategorizeDataTn(1,i) >= -200) && (CategorizeDataTn(1,i) < -150))) 
        CountDatapoints_classificationt200 = CountDatapoints_classificationt200+1;
    end
    if (((CategorizeDataTn(1,i) <= 250) && (CategorizeDataTn(1,i) > 200)) || ((CategorizeDataTn(1,i) >= -250) && (CategorizeDataTn(1,i) < -200))) 
        CountDatapoints_classificationt250 = CountDatapoints_classificationt250+1;
    end
    if (((CategorizeDataTn(1,i) <= 300) && (CategorizeDataTn(1,i) > 250)) || ((CategorizeDataTn(1,i) >= -300) && (CategorizeDataTn(1,i) < -250))) 
        CountDatapoints_classificationt300 = CountDatapoints_classificationt300+1;
    end
    if ((CategorizeDataTn(1,i) > 300) || (CategorizeDataTn(1,i) < -300))
        CountDatapoints_classificationthigher = CountDatapoints_classificationthigher+1;
    end
end
folder_id = [Maschinendaten.Entwurf.Optionen.folder_id];
file_id = [Maschinendaten.Entwurf.Optionen.file_id];

% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classificationt50,'ClassBetr','C19:C19')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classificationt100,'ClassBetr','C20:C20')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classificationt150,'ClassBetr','C21:C21')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classificationt200,'ClassBetr','C22:C22')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classificationt250,'ClassBetr','C23:C23')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classificationt300,'ClassBetr','C24:C24')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_classificationthigher,'ClassBetr','C25:C25')

%% Save results to Excel

folder_id = [Maschinendaten.Entwurf.Optionen.folder_id]; %file_id, '_data'];
file_id = [Maschinendaten.Entwurf.Optionen.file_id];

Results.E = Results.E(end)/1000/3600;
Results.S = Results.S(end);

%EConsCyc_export = struct2table(Results);
CycEff_export = array2table(meanEffMix);
% EConsCyc_export = export_excel(struct2table(Results));
% CycEff_export = export_excel(array2table(meanEffMix));
%copyfile('Ergebnisse/Vorlage_Ergebnisse.xlsx',['Ergebnisse/Ergebnisse_',file_id,'.xlsx']);
%writetable(EConsCyc_export,['3_Ergebnisse/',folder_id,'/1_Entwurf_','/Entwurf_',file_id,'.xlsx'],'Sheet','EnergyCons','Range','C3:C7','WriteVariableNames',0,'WriteRowNames',0)
%writetable(CycEff_export,['3_Ergebnisse/',folder_id,'/1_Entwurf_', '/Entwurf_',file_id,'.xlsx'],'Sheet','EnergyCons','Range','C7:C8','WriteVariableNames',0,'WriteRowNames',0)
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'],(struct2array(Results))','EnergyCons','C3:C6')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], cellstr(Cycle),'EnergyCons','C7:C7')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'],(table2array(CycEff_export)*100)','EnergyCons','C8:C8') % Umwandlung in Prozent
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_overload, 'EnergyCons','C9:C9')
% xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'], CountDatapoints_nominal,'EnergyCons','C10:C10')

FZG_LDS1 = 'Please Wait...';
disp(FZG_LDS1);



%% Start Costmodel Angerer

[RES, Optimierer] = Powertrain_Calculation(folder_id, file_id, LDS_values, Maschinendaten, vehicle, CycEff_export)

    % Start Visualization if required
    if visualize_ECS == 1
       [Tn] =  dynamicplot_motor(vehicle, driving_cycle, data_R.Tn);
    end
    

end
end
