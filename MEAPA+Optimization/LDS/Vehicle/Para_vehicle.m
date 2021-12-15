function [vehicle] = Para_vehicle(architecture, LDS_values, FZD_LDS)
% This function initializes a vehicle concept for the longitudinal dynamic
% simulation. Use this function to set all parameters according to the
% vehicle concept. 

%% Inputs
% architecture: name for a powertrain archticture. Specify the number of gearboxes and
%               motors at each axle. For example:
%   - GM_X:     one gearbox, one motor at the front axle
%   - GM_GM:    one gearbox and one motor at each axle
%   - 2G2M_X    two gearboxes, two motors at the front axle

%% Outputs
%   vehicle:    struct with all parameters of the vehicle concept 

%% Information
% Set all values according to the vehicle concept. Check the description
% after the parameters for more information on "how to obtain the values".

%% Environment
vehicle.environment.roh_L   = 1.2041;       % Density of air in kg/m^3 (sea level, 20 °C)
vehicle.environment.g       = 9.81;         % in N/kg

%% Vehicle parameters
vehicle.m               = str2double(FZD_LDS.fz_m);                  % in kg - Use Mass_calculation() regression functions 
% vehicle.maxLoad         = 500;                                     % in kg - Only important for gradeability calculation under maximum load
vehicle.rotatingMass    = vehicle.m*0.07;                            % in kg - estimation of inertia from rotating parts
vehicle.c_w             = str2double(FZD_LDS.cW);                    % in - Use regression functions in Excel sheet - english cd value
vehicle.A_front         = str2double(FZD_LDS.A);                     % in m^2 -  Use: "0.85*Width*Height" to estimate it - english reference area
vehicle.r_tyre          = str2double(FZD_LDS.tyre_r);                % in m                    
vehicle.f_R             = 0.02;                                      % in - Estimation
vehicle.battery_cap     = str2double(FZD_LDS.battcap);               % in kWh - Use battery design description
vehicle.auxiliary       = str2double(FZD_LDS.aux);                   % in W - Estimation

vehicle.Maschinendaten = LDS_values;

%% Powertrain parameters
% Use your para_powertrain_XX_XX function to initialize a powertrain with
% your specifications. Note: "crtl + d" doesn't work here, you have to open
% your function from the folder path.
Para_powertrain = str2func(strcat('Para_powertrain_', architecture)); %set Para_powertrain() to be your function 
vehicle = Para_powertrain(vehicle, FZD_LDS); 

% Instead of initializing the vehicle everytime you run the function, you
% can save it as a .mat-file after the first time and then just load it. 
function_name = '\Para_vehicle';
base_path = mfilename('fullpath');
base_path = erase(base_path,function_name);
save([base_path '\Saved_Vehicles\vehicle_eGolf'], 'vehicle'); %save the variable "vehicle" under the vehicle_name you specifiy
% load('vehicle_name');

end
