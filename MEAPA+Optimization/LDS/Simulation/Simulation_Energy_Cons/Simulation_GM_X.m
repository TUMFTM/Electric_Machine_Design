%% Simulation vehicle

function [P_EL_MIN, DATA_R]=Simulation_GM_X(V, v, a, vehicle)

%SIMULATION_GM_X calculates the required energy consumption of a vehicle
%with the powertrain GM_X (Gearbox and Motor at one axle, other axle
%powerless)

% Author:   Alexander Koch, FTM, TUM
% Date:     12/10/2018

%% Inputs
%   V:             Matrix with required Torque of vehicle in first line and revolutional speed of wheels in second line
%   vehicle:    struct with the vehicle parameters


%% Outputs
%   P_EL_MIN:      Vektor with minimal electrical power requirement for every time step
%   DATA_R:         Struct with Data about the motor torque, revolutional speed and Gear



%% Powertrain

 [V_Gearbox_Front]= Partition_Gearbox(vehicle.GEARBOX{1,1}, V);

 [P_EL_Motor_Front, DATA_Front] = Determination_Power (V_Gearbox_Front, vehicle.MOTOR{1,1}, vehicle);
 
 [P_EL_Gearbox_Front, DATA] = Determination_Gear (P_EL_Motor_Front, vehicle.GEARBOX{1,1}, DATA_Front);
 
 P_EL_MIN=P_EL_Gearbox_Front+vehicle.auxiliary; % P_EL_Gearbox stellt Maschinenleistung dar


%% Battery (Simple Eff)

%[P_EL_MIN]=BatteryEfficiency(P_EL_MIN);
 if P_EL_MIN>=0
     P_EL_MIN=P_EL_MIN/0.97;
 else
     P_EL_MIN=P_EL_MIN*0.97;
 end


%% Rearangement of DATA

% if v and a == 0 No Torque is required
DATA(1,(v==0&a==0))=0;

 DATA_R.Tn=zeros(2,length(DATA(1,:)));
 
 DATA_R.Tn(1,:)=DATA(1,:);
 DATA_R.Tn(2,:)=DATA(2,:);
 DATA_R.Gear_1=DATA(3,:);

end
 
 
