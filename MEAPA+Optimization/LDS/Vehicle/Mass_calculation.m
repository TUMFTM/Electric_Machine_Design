function [m_veh] = Mass_calculation(V_red, n_cells, P_nom, n_gears, n_pass)
%MASS_CALCULATION uses Regression functions to determine the vehicle
%weight. 

%% Inputs
% V_red:    reduced volume of the vehicle in m³. See Fuchs 2014, p. 39
% n_cells:  number of battery cells in -.
% P_nom:    nominal Power of the electric motor in kW.
% n_gears:  number of gears of the transmission in -.
% n_pass:   number of passengers in the vehicle in -.

%% Outputs
% m_veh:    estimated vehicle weight in kg

%% Information
%The functions need to be adapted before using them. Please choose a value
%for the constant parameters C1, C2, ..., Cn. 

%Important: Expand the caluclation formulas to consider multiple powertrain components.

%% Matlab tip
% Matlab allows to implement a function         f(x) = a + bx       
% as                                            f = (@x) a + bx;

%% Regression functions
m_red = @(V_red) -300 + 130 * V_red;   %reduced vehicle weight for BEVs in kg - Use these equations:
%V_red_sedan     = (0.5 * overhang_front + wheelbase + 2/3  * overhang rear) * width * height
%V_red_hatchback = (0.5 * overhang_front + wheelbase + 0.75 * overhang rear) * width * height

m_battery = @(n_cells) C3 + C4 * n_cells;                       %battery weight in kg
m_inverter = @(P_nom) C5 * P_nom;                               %inverter weight in kg
m_motor = @(P_nom) C6 * P_nom;                                  %motor weight in kg
m_transmission = @(P_nom, n_gears) C7 * P_nom + C8 * n_gears;   %transmission weight in kg
m_passenger = @(n_pass) C9 * n_pass;                            %passenger weight in kg

%% Calculation
m_veh = m_red(V_red) + m_battery(n_cells) + m_inverter(P_nom) + m_motor(P_nom) ...
    + m_transmission(P_nom, n_gears) + m_passenger(n_pass);     %estimated vehicle weight in kg
end

