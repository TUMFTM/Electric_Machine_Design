function [P_EL_MIN, E_EL_MIN, data_R, F_Ges, v, V] = Simulation_energy_consumption(dc, vehicle, architecture)

%% Calculation Resistance
v = dc.speed';                     % determine velocity vector, v in m/s
% a = [0 diff(v)];                 % determine acceleration, a in m/s²
v = (v(1:end-1)+v(2:end))/2;       % average speed between two points
a = [diff(dc.speed')];             % determine acceleration, a in m/s²
dt = diff(dc.time');

alpha = zeros(size(v));     %elevation angle, in this case 0° for all steps

[v, V, F_Ges]=Calc_Resistance(v,a,alpha, vehicle);

%% Calculation of Power requirement - different for each powertrain architecture

 Simulation=str2func(strcat('Simulation_', architecture));
 [P_EL_MIN, data_R]=Simulation(V, v, a, vehicle);
 
 E_EL_MIN=P_EL_MIN.*dt;
 
end

