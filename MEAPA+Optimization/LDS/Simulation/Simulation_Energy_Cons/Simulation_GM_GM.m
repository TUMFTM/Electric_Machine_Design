%% Simulation vehicle


function [P_EL_MIN, DATA_R]=Simulation_GM_GM(V, v, a, vehicle)

%SIMULATION_GM_GM calculates the required energy consumption of a vehicle
%with the powertrain GM_GM (Gearbox and Motor at each axle)

% Author:   Alexander Koch, FTM, TUM
% Date:     12/10/2018

%% Inputs
%   V:             Matrix with required Torque of vehicle in first line and revolutional speed of wheels in second line
%   vehicle:    struct with the vehicle parameters


%% Outputs
%   P_EL_MIN:      Vektor with minimal electrical power requirement for every time step
%   DATA_R:         Struct with Data about the motor torque, revolutional speed and Gear

vehicle.troque_range{1} = [vehicle.MOTOR{1}.eff_T_axis(end)*max(vehicle.GEARBOX{1}.gear_ratio);vehicle.MOTOR{1}.eff_T_axis(1)*min(vehicle.GEARBOX{1}.gear_ratio)];
vehicle.troque_range{2} = [vehicle.MOTOR{1}.eff_T_axis(end)*max(vehicle.GEARBOX{1}.gear_ratio);vehicle.MOTOR{1}.eff_T_axis(1)*min(vehicle.GEARBOX{1}.gear_ratio)];


%% Powertrain

 [V_Axle_Front, V_Axle_Rear]= Partition_Torque(vehicle.troque_range{1}, vehicle.troque_range{2}, V);
 
 [V_Gearbox_Front]= Partition_Gearbox(vehicle.GEARBOX{1,1}, V_Axle_Front);
 
 [V_Gearbox_Rear]= Partition_Gearbox(vehicle.GEARBOX{1,1},V_Axle_Rear);

 [P_EL_Motor_Front, DATA_Front] = Determination_Power (V_Gearbox_Front, vehicle.MOTOR{1,1});

 [P_EL_Motor_Rear, DATA_Rear] = Determination_Power (V_Gearbox_Rear, vehicle.MOTOR{1,1});
 
 [P_EL_Getriebe_Front, DATA_Front] = Determination_Gear (P_EL_Motor_Front, vehicle.GEARBOX{1,1}, DATA_Front);

 [P_EL_Getriebe_Rear, DATA_Rear] = Determination_Gear (P_EL_Motor_Rear, vehicle.GEARBOX{1,1}, DATA_Rear); 
 
 [P_EL_MIN, DATA] = Determination_Distribution (P_EL_Getriebe_Front,P_EL_Getriebe_Rear, [DATA_Front;DATA_Rear]);
 
 P_EL_MIN=P_EL_MIN+vehicle.auxiliary;
  

%% Battery (Simple Eff)
[P_EL_MIN]=BatteryEfficiency(P_EL_MIN);

%% Rearangement of DATA

% if v and a == 0 No Torque is required
DATA(1,(v==0&a==0))=0;
DATA(4,(v==0&a==0))=0;


 DATA_R.Tn=zeros(4,length(DATA(1,:)));
 
 
 DATA_R.Tn(1,:)=DATA(1,:);
 DATA_R.Tn(2,:)=DATA(2,:);
 DATA_R.Tn(3,:)=DATA(4,:);
 DATA_R.Tn(4,:)=DATA(5,:);
 DATA_R.Gear_1=DATA(3,:);
 DATA_R.Gear_2=DATA(6,:);
 
 
end
 
 
