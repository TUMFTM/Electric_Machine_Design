function [vehicle] = Para_powertrain_GM_GM(vehicle, FZD_LDS)
%PARA_POWERTRAIN_GM_X initializes a powertrain configuration. Adapt this
%function for different powertrain architectures.

%% Matlab Tip:
% In this example, only one front motor and gearbox is parameterized.
% However, the code for the other elements is already implemented. In
% Matlab, you set lines of code to "comment" using "crtl + R". Use 
% "ctrl + T" to uncomment. Afterwards, just adapt the given values
% according to your powertrain design.

%% Powertrain
%initialize - gearbox and motor always need to have 4 elements. If a
%powertrain contains less motors (and gearboxes), the fields remain
%empty.
vehicle.GEARBOX = cell(1,4); %initialize empty gearbox struct
vehicle.MOTOR = cell(1,4); %initialize empty motor struct

%% -----------------------    Front 1    ----------------------- 
% Gearbox front 1
% gear ratios of each gear in -; [1.gear 2.gear ... n.gear], Use transmission design guide to determine the ratios
vehicle.GEARBOX{1}.gear_ratio       = [str2double(FZD_LDS.GearRatio)];         
vehicle.GEARBOX{1}.eff              = [0.96]; % efficiencies in -;  [eff1.gear  eff2.gear ... effn.gear], besser 0,97 laut [CHO13]

motor=vehicle.Maschinendaten; %Übergabe

% Motor front 1
motor.type       = vehicle.Maschinendaten.type;          % PSM or ASM
motor.power      = vehicle.Maschinendaten.power;         % in kW
motor.n_n        = vehicle.Maschinendaten.n_n;           % in 1/min
motor.n_max      = vehicle.Maschinendaten.n_max;         % in 1/min
motor.U_n        = vehicle.Maschinendaten.U_n;           % in V
% Generate the motor data with motor tool:
vehicle.MOTOR{1}=vehicle.Maschinendaten;
% vehicle.MOTOR{1} = GENERATE_motor_efficiency_map(motor); %Generiert Map aus Horlbecktool, wenn bereits Motor vorhanden, nicht notwendig


% %% -----------------------    Front 2    ----------------------- 
% % Gearbox front 2
% % gear ratios of each gear in -; [1.gear 2.gear ... n.gear], Use transmission design guide to determine the ratios
% vehicle.GEARBOX{2}.gear_ratio       = ;        
% vehicle.GEARBOX{2}.eff              = ; % efficiencies in -;  [eff1.gear  eff2.gear ... effn.gear]
% 
% % Motor front 2
% motor.type       = ;        % PSM or ASM
% motor.power      = ;           % in kW
% motor.n_n        = ;         % in 1/min
% motor.n_max      = ;        % in 1/min
% motor.U_n        = ;          % in V
% % Generate the motor data with motor tool:
% vehicle.MOTOR{2} = GENERATE_motor_efficiency_map(motor);


%% -----------------------    Rear 1    ----------------------- 
% Gearbox rear 1
% gear ratios of each gear in -; [1.gear 2.gear ... n.gear], Use transmission design guide to determine the ratios
vehicle.GEARBOX{1}.gear_ratio       = [str2double(FZD_LDS.GearRatio)];         
vehicle.GEARBOX{1}.eff              = [0.96]; % efficiencies in -;  [eff1.gear  eff2.gear ... effn.gear], besser 0,97 laut [CHO13]

motor=vehicle.Maschinendaten; %Übergabe

% Motor rear 1
motor.type       = vehicle.Maschinendaten.type;          % PSM or ASM
motor.power      = vehicle.Maschinendaten.power;         % in kW
motor.n_n        = vehicle.Maschinendaten.n_n;           % in 1/min
motor.n_max      = vehicle.Maschinendaten.n_max;         % in 1/min
motor.U_n        = vehicle.Maschinendaten.U_n;           % in V
% Generate the motor data with motor tool:
vehicle.MOTOR{1}=vehicle.Maschinendaten;
% vehicle.MOTOR{1} = GENERATE_motor_efficiency_map(motor); %Generiert Map aus Horlbecktool, wenn bereits Motor vorhanden, nicht notwendig


% %% -----------------------    Rear 2    ----------------------- 
% % Gearbox rear 2
% % gear ratios of each gear in -; [1.gear 2.gear ... n.gear], Use transmission design guide to determine the ratios
% vehicle.GEARBOX{4}.gear_ratio       = ; 
% vehicle.GEARBOX{4}.eff              = ; % efficiencies in -;  [eff1.gear  eff2.gear ... effn.gear]
% 
% % Motor rear 2
% motor.type       = ;        % PSM or ASM
% motor.power      = ;           % in kW
% motor.n_n        = ;         % in 1/min
% motor.n_max      = ;        % in 1/min
% motor.U_n        = ;          % in V
% % Generate the motor data with motor tool:
% vehicle.MOTOR{4} = GENERATE_motor_efficiency_map(motor);


end

