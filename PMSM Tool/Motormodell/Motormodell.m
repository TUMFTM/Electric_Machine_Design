% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Model for the design of a permanent magnet excited synchronous machine and
% subsequent efficiency map calculation
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (kalt@ftm.mw.tum.de)
%           Jonathan Erhard
%           Prof. Markus Lienkamp
% -------------------------------------------------------------------------

function [i_d_mesh, i_q_mesh, M_max_vec, M_max_mesh, omega_k_mesh] = Motormodell(prim, ent, reg, omega_k_vec, tics_M)
% Motor model - This function calculates the full load characteristic and the currents
% and voltages of a PMSM.
% In the motor model function, the full load characteristic curve is first displayed under
% observance of current and voltage limits. This is then used to calculate
% the currents and voltages for operation in the permissible range

%% Create torque vector
M_vec = linspace(0, ent.M_N, tics_M);

%% Create speed/torque grid
[omega_k_mesh, M_mesh] = meshgrid(omega_k_vec, M_vec);

%% Calculation of full load curve (calculate maximum torque for speed vector)
% Maximum torque under auxiliary conditions
[i_s, beta] = Optimierung_M(prim, ent, reg, omega_k_mesh);
M_max_vec = 1.5.*prim.p.*((ent.L_d-ent.L_q).*i_s.^2.*0.5.*sin(2.*beta) + ent.psi_PM.*i_s.*sin(beta));

%% Create torque vector
M_vec = linspace(0, max(M_max_vec), tics_M);

%% Create speed/torque grid
[omega_k_mesh, M_max_mesh] = meshgrid(omega_k_vec, M_vec);

%% Adjustment of map range to full load curve
% Tolerance for max torque (numerical inaccuracy)
tol = 1e-4;
% Filtering operating points that cannot be reached
for i = 1:length(omega_k_mesh(1,:))
    M_max_mesh(M_max_mesh(:,i)>M_max_vec(i)+tol,i) = NaN;
end

%% Calculation of currents for given speeds and torques
% Minimization of currents under auxiliary conditions
[i_d_mesh, i_q_mesh] = Optimierung_i(prim, ent, reg, omega_k_mesh, M_max_mesh);