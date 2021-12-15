function [RESacc, v_acc] = Simulation_acceleration(vehicle, v_max_sim, t_sim)
% SIMULATION_ACCELERATION determines the maximum acceleration of the vehicle concept 

% Author:   Sebastian Krapf, FTM, TUM
% Date:     02/10/2018

%% Inputs
%   vehicle:    struct with the vehicle parameters
%   v_max_sim:  max. simulated velocity in km/h. If vehicle reaches this velocity,
%               the simulation ends
%   t_sim:      simulation time in s. Simulation ends after this time
%               automatically

%% Outputs
%   v_acc:      struct with a vehicle speed vector in km/h 
%               and a time vector in s

tic
%% Initialization
%siulation parameters
t_max = t_sim;      %maximum simulation time in s
delta_t = 0.01;     %simulation step in s
%velocity vector
v = 0:1:v_max_sim;      %vehicle velocity in km/h;

%% Calculation
%% Maximum possible torque of each motor for each velocity point (including gearbox efficiency)
T_max_wheel_v_all = zeros(size(v));
for i = 1:numel(vehicle.MOTOR)
    if ~isempty(vehicle.GEARBOX{i})
    % calculate the motor speed for each gear at each vehicle speed 
    n_wheel = v/3.6/vehicle.r_tyre *60/2/pi;     %in rpm, rotational speed wheel
    n_motor_v = (n_wheel'*vehicle.GEARBOX{i}.gear_ratio);%'; %in rpm, rotational speed motor for each gear
    
    % determine the respective maximum torque --> T_max(motor,v,i)
    for j = 1:length(n_motor_v,1)
        F = griddedInterpolant(vehicle.MOTOR{i}.T_max_n_axis,vehicle.MOTOR{i}.T_max,'linear','none'); %interpolate values from the T_max motor curve
        T_motor_v(j,:) = F(n_motor_v(j,:)); %use griddedInterpolant
    end
    T_wheel_v = T_motor_v.*vehicle.GEARBOX{i}.gear_ratio'.*vehicle.GEARBOX{i}.eff'; %maximum torque at wheel depending on the gear ratio
    
    % calculate maximum torque at the wheel depending on vehicle speed and gear to reach the maximum torque
    [T_max_wheel_v(i,:), id_gear] = max(T_wheel_v,[],1);
    end
end

%% Sum of all possible torques at each velocity point
T_max_wheel_v_all =sum(T_max_wheel_v,1) ; %vector with the max. torque of all drives

%% Forward Simulation: Calculate resulting acceleration with maximum torque 
t = 1;
v_actual = 0;
F = griddedInterpolant(v,T_max_wheel_v_all,'linear','none'); %interpolation function

while (t < t_max/delta_t && v_actual(t) < max(v))
       

    % maximum available torque at wheels from all motors
    T = F(v_actual(t)); %use interpolation function to interpolate
    
    % resulting acceleration with this torque (fundamental equation of vehicle longitudinal dynamics)
    F_drive = T / vehicle.r_tyre;
    F_aero = 0.5 * vehicle.environment.roh_L * vehicle.c_w * vehicle.A_front * (v_actual(t)/3.6).^2; 
    F_fric = vehicle.m * vehicle.environment.g * vehicle.f_R * cos(0); %no elevation (alpha = 0)
    F_slope = 0;
    F_acc = F_drive -  F_aero - F_fric - F_slope;
    
    acc(t) = F_acc/(vehicle.m+vehicle.rotatingMass); %vector with the acceleration of each simulation step
    
    % resulting velocity of next time step
    v_actual(t+1) = v_actual(t)+ acc(t)*3.6 * delta_t; %vector with the veloctiy of each simulation step
    t = t+1; % count timer up
    
        % break if v_actual = NaN: This means n_mot > n_mot_max
 if isnan(v_actual(t))
        disp('Vehicle reached maximum motor rotational speed. End of simulation.');
        break;
        end
end

v_acc.v = v_actual;                 %struct that stores result: velocity vector
v_acc.t = 0:delta_t:delta_t*(t-1);  %struct that stores result: time vector

RESacc = struct;
RESacc.v_acc = v_acc;
RESacc.F_drive = F_drive;
RESacc.F_aero = F_aero;
RESacc.F_fric = F_fric;
RESacc.F_slope = F_slope;
RESacc.F_acc = F_acc;

toc

end
