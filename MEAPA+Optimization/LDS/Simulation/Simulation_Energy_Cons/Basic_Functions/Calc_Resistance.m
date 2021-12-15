
function [v, V, F_Ges]=Calc_Resistance(v,a,alpha, vehicle)


%% Calculation of Resistance

F_Air=0.5*vehicle.environment.roh_L*vehicle.c_w*vehicle.A_front*v.^2;
F_Roll=vehicle.m*vehicle.environment.g*vehicle.f_R*cos(alpha);
F_Slope=vehicle.m*vehicle.environment.g*sin(alpha);
F_Acc=(vehicle.m+vehicle.rotatingMass)*a;

F_Ges= F_Air + F_Roll + F_Slope + F_Acc;


%% Calculation of required Torque and revs
Torque_Total=F_Ges*vehicle.r_tyre;

n= v/vehicle.r_tyre/2/pi;   % 1/s

V=[Torque_Total;n];

% Antriebsleistung P = F_Ges *v
end