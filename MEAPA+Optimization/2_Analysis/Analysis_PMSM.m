% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Modell for the design and analysis of PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (svenja.kalt@tum.de), 
%           Jonathan Erhard 
% -------------------------------------------------------------------------


%% Table of contents
% A) Preprocessing
% B) Torque control
%   B.1) Motormodel
%   B.2) Generatormodel
% C) Loss calculation
% D) Map creation
% E) Postprocessing
% F) Help functions
%   F.1) Optimierung_M %Optimization Torque
%   F.2) Optimierung_i %Optimization Current

function [Analyse] = Analysis_PMSM(handles)
%Analysis This function performs the analysis of a PMSM.
% In the function Analysis_PMSM maps are calculated for the previously
% designed machine. First the currents and voltages are
% voltages are calculated by means of a motor model. Then the
% individual loss components are determined. Finally the
% efficiency map is generated. Note: There is no thermal
% check or similar. The user is therefore advised to feed the machine only with realistic
% with realistic data (i_max and u_max).
% V/A: linear & stationary consideration:
% - no saturation
% - no cross coupling of inductors
% - inductors are not current dependent
% - no dynamics
% - chained flow of PM constant

%% A) Preprocessing
% #########################################################################
% #   A) PREPROCESSING                                                    #
% #########################################################################

% Restoring the variables for easier use
const = handles.Entwurf.Konstanten;
emag = handles.Entwurf.EMAG;
geo = handles.Entwurf.Geometrie;
wick = handles.Entwurf.Wicklung;

rated = handles.Entwurf.Bemessungswerte;
richt = handles.Entwurf.Richtwerte;
opt_Entwurf = handles.Entwurf.Optionen;

opt = handles.opt;

%% B) Torque control
% #########################################################################
% #   B) TORQUE CONTROL                                                   #
% #########################################################################

% Vector rotaional speed
betr.omega_k_max = (opt.n_max/60)*(2*pi)*rated.p;

i_f = emag.Psi_PM / emag.L_d;
if(i_f>=opt.i_1max)
    omega_k_max2 = sqrt((opt.u_1max^2 - emag.R_1^2*opt.i_1max^2) / (emag.Psi_PM - emag.L_d*opt.i_1max)^2) - 4;
    if(betr.omega_k_max>omega_k_max2)
        betr.omega_k_max = omega_k_max2;
        opt.n_max = betr.omega_k_max ./ (rated.p .* (2.*pi)) .* 60;
        warning('max. speed limited by high interlinked flow of PM');
    end
end
clear i_f omega_k_max2

betr.omega_k_vec = linspace(0, betr.omega_k_max, opt.n_tics);

%% B.1) Motormodel
% #########################################################################
% #   B.1) MOTORMODEL                                                     #
% #########################################################################

% Calculation of the full load characteristic (calculate maximum torque for speed vector).
% Maximization of the torque under secondary conditions (operating range motor)
[ctrl.mot.i_1, ctrl.mot.beta] = Optimierung_M(rated, emag, opt, betr.omega_k_vec, 0);
betr.mot.M_max_vec = 1.5.*rated.p.*((emag.L_d-emag.L_q).*ctrl.mot.i_1.^2.*0.5.*sin(2.*ctrl.mot.beta) + emag.Psi_PM.*ctrl.mot.i_1.*sin(ctrl.mot.beta));
betr.mot.M_max_vec = betr.mot.M_max_vec';

% Tolerance for max. moment (reason: numerical inaccuracy)
tol = 1e-3;

% Create torque vector
betr.mot.M_vec = linspace(max(betr.mot.M_max_vec)-tol, 0, opt.M_tics);
clear tol

% Create speed/torque grid
[betr.mot.omega_k_mesh, betr.mot.M_max_mesh] = meshgrid(betr.omega_k_vec, betr.mot.M_vec);

% Adaptation of characteristic map range to full load characteristic curve
% Filtering of the points that cannot be approached
for i = 1:length(betr.mot.omega_k_mesh(1,:))
    betr.mot.M_max_mesh(betr.mot.M_max_mesh(:,i)>betr.mot.M_max_vec(i),i) = NaN;
end


% Calculation of currents (stator) for points in the speed/torque grid.
% Minimization of currents under constraints (motor operating range)
[ctrl.mot.i_1d_mesh, ctrl.mot.i_1q_mesh] = Optimierung_i(rated, emag, opt, betr.mot.omega_k_mesh, betr.mot.M_max_mesh, 0);

% Calculation of the voltages (stator) for points in the speed/torque grid.
ctrl.mot.u_1d_mesh = emag.R_1.*ctrl.mot.i_1d_mesh - betr.mot.omega_k_mesh.*emag.L_q.*ctrl.mot.i_1q_mesh;
ctrl.mot.u_1q_mesh = emag.R_1.*ctrl.mot.i_1q_mesh + betr.mot.omega_k_mesh.*(emag.L_d.*ctrl.mot.i_1d_mesh + emag.Psi_PM);

% Calculation of electric power
ctrl.mot.P_el_mesh = 1.5.*(ctrl.mot.i_1d_mesh.*ctrl.mot.u_1d_mesh + ctrl.mot.i_1q_mesh.*ctrl.mot.u_1q_mesh);

% Calculation speed vectors
betr.mot.omega_el_mesh = betr.mot.omega_k_mesh; 
betr.mot.f_el_mesh = betr.mot.omega_el_mesh ./ (2.*pi);  
betr.mot.omega_m_mesh = (betr.mot.omega_el_mesh ./ rated.p);
betr.mot.n_m_mesh = (betr.mot.omega_m_mesh ./ (2.*pi)) .* 60; 
betr.mot.n_m_vec = betr.mot.n_m_mesh(1,:);

%% B.2) Generatormodel
% #########################################################################
% #   B.2) GENERATORMODEL                                                 #
% #########################################################################

if(opt.Generator)
    % Calculation of the full load characteristic (calculate maximum torque for speed vector)
    % Maximization of the torque under secondary conditions (operating range generator)
    [ctrl.gen.i_1, ctrl.gen.beta] = Optimierung_M(rated, emag, opt, betr.omega_k_vec, 1);
    betr.gen.M_max_vec = 1.5.*rated.p.*((emag.L_d-emag.L_q).*ctrl.gen.i_1.^2.*0.5.*sin(2.*ctrl.gen.beta) + emag.Psi_PM.*ctrl.gen.i_1.*sin(ctrl.gen.beta));
    betr.gen.M_max_vec = betr.gen.M_max_vec';
    
    % Tolerance for max moment (reason: numerical inaccuracy)
    tol = 1e-3;
    
    % Create torque vector
    betr.gen.M_vec = linspace(0, min(betr.gen.M_max_vec)+tol, opt.M_tics);
    clear tol
    
    % Create speed/torque grid
    [betr.gen.omega_k_mesh, betr.gen.M_max_mesh] = meshgrid(betr.omega_k_vec, betr.gen.M_vec);
    
    % Adaptation of characteristic map range to full load characteristic curve
    % Filtering of the points that cannot be approached
    for i = 1:length(betr.gen.omega_k_mesh(1,:))
        betr.gen.M_max_mesh(betr.gen.M_max_mesh(:,i)<betr.gen.M_max_vec(i),i) = NaN;
    end
    
    % Calculation of currents (stator) for points in speed/torque grid
    % Minimization of currents under constraints (generator operating range)
    [ctrl.gen.i_1d_mesh, ctrl.gen.i_1q_mesh] = Optimierung_i(rated, emag, opt, betr.gen.omega_k_mesh, betr.gen.M_max_mesh, 1);

    % Calculation voltages (stator) for points in speed/torque grid
    ctrl.gen.u_1d_mesh = emag.R_1.*ctrl.gen.i_1d_mesh - betr.gen.omega_k_mesh.*emag.L_q.*ctrl.gen.i_1q_mesh;
    ctrl.gen.u_1q_mesh = emag.R_1.*ctrl.gen.i_1q_mesh + betr.gen.omega_k_mesh.*(emag.L_d.*ctrl.gen.i_1d_mesh + emag.Psi_PM);
    
    % Calculation of the el. power
    ctrl.gen.P_el_mesh = 1.5.*(ctrl.gen.i_1d_mesh.*ctrl.gen.u_1d_mesh + ctrl.gen.i_1q_mesh.*ctrl.gen.u_1q_mesh);

    % Calculation speed vectors
    betr.gen.omega_el_mesh = betr.gen.omega_k_mesh; % elektrische Winkelgeschwindigkeit
    betr.gen.f_el_mesh = betr.gen.omega_el_mesh ./ (2.*pi);  % elektrische Frequenz
    betr.gen.omega_m_mesh = (betr.gen.omega_el_mesh ./ rated.p); % mechanische Winkelgeschwindigkeit
    betr.gen.n_m_mesh = (betr.gen.omega_m_mesh ./ (2.*pi)) .* 60; % mechanische Drehzahl
    betr.gen.n_m_vec = betr.gen.n_m_mesh(1,:);
    
    % Combine engine and generator
    % ctrl.i_1 = [ctrl.mot.i_1; ctrl.gen.i_1(2:end,:)];
    % ctrl.beta = [ctrl.mot.beta; ctrl.gen.beta(2:end,:)];
    % betr.M_max_vec = [betr.mot.M_max_vec betr.gen.M_max_vec(:,2:end)];
    betr.M_vec = [betr.mot.M_vec betr.gen.M_vec(:,2:end)];
    betr.omega_k_mesh = [betr.mot.omega_k_mesh; betr.gen.omega_k_mesh(2:end,:)];
    betr.M_max_mesh = [betr.mot.M_max_mesh; betr.gen.M_max_mesh(2:end,:)];
    ctrl.i_1d_mesh = [ctrl.mot.i_1d_mesh; ctrl.gen.i_1d_mesh(2:end,:)];
    ctrl.i_1q_mesh = [ctrl.mot.i_1q_mesh; ctrl.gen.i_1q_mesh(2:end,:)];
    ctrl.u_1d_mesh = [ctrl.mot.u_1d_mesh; ctrl.gen.u_1d_mesh(2:end,:)];
    ctrl.u_1q_mesh = [ctrl.mot.u_1q_mesh; ctrl.gen.u_1q_mesh(2:end,:)];
    ctrl.P_el_mesh = [ctrl.mot.P_el_mesh; ctrl.gen.P_el_mesh(2:end,:)];
    betr.omega_el_mesh = [betr.mot.omega_el_mesh; betr.gen.omega_el_mesh(2:end,:)];
    betr.f_el_mesh = [betr.mot.f_el_mesh; betr.gen.f_el_mesh(2:end,:)];
    betr.omega_m_mesh = [betr.mot.omega_m_mesh; betr.gen.omega_m_mesh(2:end,:)];
    betr.n_m_mesh = [betr.mot.n_m_mesh; betr.gen.n_m_mesh(2:end,:)];
    % betr.n_m_vec = [betr.mot.n_m_vec betr.gen.n_m_vec(:,2:end)];
else
    % Combine engine and generator
    % ctrl.i_1 = ctrl.mot.i_1;
    % ctrl.beta = ctrl.mot.beta;
    % betr.M_max_vec = betr.mot.M_max_vec;
    betr.M_vec = betr.mot.M_vec;
    betr.omega_k_mesh = betr.mot.omega_k_mesh;
    betr.M_max_mesh = betr.mot.M_max_mesh;
    ctrl.i_1d_mesh = ctrl.mot.i_1d_mesh;
    ctrl.i_1q_mesh = ctrl.mot.i_1q_mesh;
    ctrl.u_1d_mesh = ctrl.mot.u_1d_mesh;
    ctrl.u_1q_mesh = ctrl.mot.u_1q_mesh;
    ctrl.P_el_mesh = ctrl.mot.P_el_mesh;
    betr.omega_el_mesh = betr.mot.omega_el_mesh;
    betr.f_el_mesh = betr.mot.f_el_mesh;
    betr.omega_m_mesh = betr.mot.omega_m_mesh;
    betr.n_m_mesh = betr.mot.n_m_mesh;
    % betr.n_m_vec = betr.mot.n_m_vec;

end

%% C) Loss calculation
% #########################################################################
% #   C) LOSS CALCULATION                                                 #
% #########################################################################

% Calculation of winding losses P_vw [W] (copper losses/current heat losses)
% Source: [Mueller08, S.438 - Formel 6.3.18]
% the 1.5 come from the Clarke transformation
% (i_d_mesh.^2 + i_q_mesh.^2) are peak values
if(opt.P_vw==1)
    loss.P_vw_mesh = 1.5 * emag.R_1 * (ctrl.i_1d_mesh.^2 + ctrl.i_1q_mesh.^2);
else
    loss.P_vw_mesh = zeros(size(ctrl.i_1d_mesh));
end

% Calculation of the core magnetization losses P_vu [W].

if(opt.P_vu==1)
    if(strcmp(opt.P_vu_Modell,'Modellansatz Jordan'))
        % Variant 1: Calculation of iron losses via inductions and
        % Electrical sheet data sheets
        
        % Calculate sheet parameters for B=1T and f=50Hz (stator)
        % Source: [Neuschl, p.32-33]
        loss.vu.w_1Fe = opt_Entwurf.Stator_Eisenmaterial.p_vFe(opt_Entwurf.Stator_Eisenmaterial.p_vFe_B_vec==1,:) ./ opt_Entwurf.Stator_Eisenmaterial.p_vFe_f_vec;
        [loss.vu.fitresult_1, ~] = fit(opt_Entwurf.Stator_Eisenmaterial.p_vFe_f_vec', loss.vu.w_1Fe', 'poly1');
        loss.vu.fitresult_data_1 = coeffvalues(loss.vu.fitresult_1);
        loss.vu.sigma_1h_s = loss.vu.fitresult_data_1(2);
        loss.vu.sigma_1w_s = loss.vu.fitresult_data_1(1);
        loss.vu.sigma_1h = loss.vu.sigma_1h_s * 50;
        loss.vu.sigma_1w = loss.vu.sigma_1w_s * 50^2;
        
        % Calculate sheet parameters for B=1T and f=50Hz (rotor)
        % Source: [Neuschl, p.32-33]
        loss.vu.w_2Fe = opt_Entwurf.Stator_Eisenmaterial.p_vFe(opt_Entwurf.Rotor_Eisenmaterial.p_vFe_B_vec==1,:) ./ opt_Entwurf.Stator_Eisenmaterial.p_vFe_f_vec;
        [loss.vu.fitresult_2, ~] = fit(opt_Entwurf.Stator_Eisenmaterial.p_vFe_f_vec', loss.vu.w_2Fe', 'poly1');
        loss.vu.fitresult_data_2 = coeffvalues(loss.vu.fitresult_2);
        loss.vu.sigma_2h_s = loss.vu.fitresult_data_2(2);
        loss.vu.sigma_2w_s = loss.vu.fitresult_data_2(1);
        loss.vu.sigma_2h = loss.vu.sigma_2h_s * 50;
        loss.vu.sigma_2w = loss.vu.sigma_2w_s * 50^2;
        
        % Surcharge factors
        % Source: [Neuschl, p.29 - formula 4.10]
        loss.vu.c_h = 1;
        loss.vu.c_w = 1;
        
        % Calculate specific losses in the different sections.
        % Source: [Neuschl, p.28 - Formula 4.9]
        loss.vu.p_1r = (loss.vu.sigma_1h.*loss.vu.c_h.*(betr.f_el_mesh./50) + loss.vu.sigma_1w.*loss.vu.c_w.*(betr.f_el_mesh./50).^2) .* (emag.B_1r./1).^2;
        loss.vu.p_1z = (loss.vu.sigma_1h.*loss.vu.c_h.*(betr.f_el_mesh./50) + loss.vu.sigma_1w.*loss.vu.c_w.*(betr.f_el_mesh./50).^2) .* (emag.B_1z_m./1).^2;
        loss.vu.p_2r = (loss.vu.sigma_2h.*loss.vu.c_h.*(betr.f_el_mesh./50) + loss.vu.sigma_2w.*loss.vu.c_w.*(betr.f_el_mesh./50).^2) .* (emag.B_2r./1).^2;
        
        % Surcharge factors
        % Source: [Mueller08, p.452 - Table 6.4.2]
        loss.vu.k_u_r = 1.7; % between 1.5 and 1.8
        loss.vu.k_u_z = 1.8; % between 1.7 and 2.5

        % Calculate losses from the specific losses
        loss.vu.P_1r = loss.vu.k_u_r .* loss.vu.p_1r .* opt_Entwurf.Stator_Eisenmaterial.rho_Fe .* geo.Vo_1r;
        loss.vu.P_1z = loss.vu.k_u_z .* loss.vu.p_1z .* opt_Entwurf.Stator_Eisenmaterial.rho_Fe .* geo.Vo_1z;
        loss.vu.P_2r = loss.vu.k_u_r .* loss.vu.p_2r .* opt_Entwurf.Rotor_Eisenmaterial.rho_Fe .* geo.Vo_2r;
        
        % Total iron losses
        loss.P_vu_mesh = loss.vu.P_1r + loss.vu.P_1z + loss.vu.P_2r;
    elseif(strcmp(opt.P_vu_Modell,'Abschaetzung ueber verketteten Fluss'))
        % Variant 2: Calculation of iron losses via interlinked flow
        % Source: [Script Drive Control for Electric Vehicles, TUM EAL]
        % Note: not recommended, because of empirical factors very inaccurate
        
        % Empirical factors
        loss.vu.c_Fe1 = 1.5; % between 1.5 and 1.6
        loss.vu.c_Fe2 = 1.5; % between 1.5 and 1.6
        
        % Calculation of iron losses
        loss.P_vu_mesh = loss.vu.c_Fe1.*betr.omega_k_mesh.^(loss.vu.c_Fe2).*((emag.L_d.*ctrl.i_1d_mesh + emag.Psi_PM).^2 + (emag.L_q.*emag.i_1q_mesh).^2);
    else
        error('Incorrect input for variable "opt.P_vu_Model"');
    end
else
    loss.P_vu_mesh = zeros(size(ctrl.i_1d_mesh));
end

% Calculation of mechanical losses P_vme [W].
% Source: [Mueller08, p.433 - formula 6.2.1]
if(opt.P_vme==1)
    % Factors of gas and bearing friction k_rb [Ws^2/m^4]
    loss.vme.k_rb = 10; % should be between 5 and 15

    % Rotor circumferential speed [m/s]
    loss.vme.v_2 = (geo.D_2a ./ 2) .* betr.omega_m_mesh;

    loss.P_vme_mesh = loss.vme.k_rb .* geo.D_2a .* (geo.l_i + 0.8^3 .* 0.6 .* geo.tau_1p) .* loss.vme.v_2.^2;
else
	loss.P_vme_mesh = zeros(size(ctrl.i_1d_mesh));
end

% Calculation of additional losses P_zus [W].
% Source: [DIN EN 60034-2-1]
if(opt.P_vzus==1)
    loss.P_vzus_mesh = abs(ctrl.P_el_mesh) .* (0.025 - 0.005*log(rated.P_N*1e-3));
else
    loss.P_vzus_mesh = zeros(size(ctrl.i_1d_mesh));
end

% Calculation of total losses P_ges [W]
loss.P_vges_mesh = loss.P_vw_mesh + loss.P_vu_mesh + loss.P_vme_mesh + loss.P_vzus_mesh;

%% D) Map Creation
% #########################################################################
% #   D) MAP CREATION                                                     #
% #########################################################################

% Calculate the mech. power
betr.mot.P_mech_mesh = betr.mot.M_max_mesh .* betr.mot.omega_m_mesh;

% Calculation of efficiency maps (motor operating range)
[var1, ~] = size(betr.mot.P_mech_mesh);
eta.mot.eta_vw_mesh_alt = betr.mot.P_mech_mesh ./ ctrl.mot.P_el_mesh;
eta.mot.eta_vw_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vw_mesh(1:var1,1:end));
eta.mot.eta_fe_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vu_mesh(1:var1,1:end));
eta.mot.eta_vme_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vme_mesh(1:var1,1:end));
eta.mot.eta_zus_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vzus_mesh(1:var1,1:end));
eta.mot.eta_ges_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vges_mesh(1:var1,1:end));

if(opt.Generator)
    % Calculate the mech. power
    betr.gen.P_mech_mesh = betr.gen.M_max_mesh .* betr.gen.omega_m_mesh;
    
    % Calculation of efficiency maps (generator operating range)
    eta.gen.eta_vw_mesh_alt = ctrl.gen.P_el_mesh ./ betr.gen.P_mech_mesh;
    eta.gen.eta_vw_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vw_mesh(var1:end,1:end));
    eta.gen.eta_fe_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vu_mesh(var1:end,1:end));
    eta.gen.eta_vme_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vme_mesh(var1:end,1:end));
    eta.gen.eta_zus_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vzus_mesh(var1:end,1:end));
    eta.gen.eta_ges_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vges_mesh(var1:end,1:end));
    
    % Combine engine and generator
    eta.eta_vw_mesh = [eta.mot.eta_vw_mesh; eta.gen.eta_vw_mesh(2:end,:)];
    eta.eta_fe_mesh = [eta.mot.eta_fe_mesh; eta.gen.eta_fe_mesh(2:end,:)];
    eta.eta_vme_mesh = [eta.mot.eta_vme_mesh; eta.gen.eta_vme_mesh(2:end,:)];
    eta.eta_zus_mesh = [eta.mot.eta_zus_mesh; eta.gen.eta_zus_mesh(2:end,:)];
    eta.eta_ges_mesh = [eta.mot.eta_ges_mesh; eta.gen.eta_ges_mesh(2:end,:)];
else
    % Combine engine and generator
    eta.eta_vw_mesh = eta.mot.eta_vw_mesh;
    eta.eta_fe_mesh = eta.mot.eta_fe_mesh;
    eta.eta_vme_mesh = eta.mot.eta_vme_mesh;
    eta.eta_zus_mesh = eta.mot.eta_zus_mesh;
    eta.eta_ges_mesh = eta.mot.eta_ges_mesh;
end

%% E) Postprocessing
% #########################################################################
% #   E) POSTPROCESSING                                                   #
% #########################################################################

% Calculation successfull
opt.Locked = 1;

% Sorting the structs
ctrl = orderfields(ctrl);
ctrl.mot = orderfields(ctrl.mot);
if(opt.Generator)
    ctrl.gen = orderfields(ctrl.gen);
    betr.gen = orderfields(betr.gen);
end
loss = orderfields(loss);
if(isfield(loss,'vu'))
    loss.vu = orderfields(loss.vu);
end
if(isfield(loss,'vme'))
    loss.vme = orderfields(loss.vme);
end
betr = orderfields(betr);
betr.mot = orderfields(betr.mot);
eta = orderfields(eta);

opt = orderfields(opt);

% Save analysis
Analyse.Momentensteuerung = ctrl;
Analyse.Verluste = loss;
Analyse.Betriebsdaten = betr;
Analyse.Wirkungsgrad = eta;
Analyse.Optionen = opt;
Analyse.EMAG = emag;
Analyse.Geometrie = geo;
Analyse.Optionen = opt;

% Save struct analysis
%save(['3_Results/',handles.Entwurf.Optionen.folder_id,'/2_Analyse/Analyse_',handles.Entwurf.Optionen.file_id,'.mat'],'Analyse');

end

%% F) Help function
% #########################################################################
% #   F) HELP FUNCTION                                                    #
% #########################################################################

%% F.1) Optimierung_M - optimization of torque
% #########################################################################
% #   F.1) OPTIMIERUNG_M - OPTIMIZATION OF TORQUE                         #
% #########################################################################

function [i_1, beta] = Optimierung_M(rated, emag, opt, omega_vec, Generator)
% Optimization Maximizes the moment under constraints and limits.
% Notation: i_1 = x(1), beta = x(2)
% Bound (motor): 0<=i_1<=i_max, 0<=beta<=pi
% Linear Inequality Constraint: none
% Linear Equality Constraint: none
% Nonlinear Constraints: see function nonlcon_Optimierung_M

if(Generator)
    % Minimization function 
    fun = @(x) -1.5.*rated.p.*((emag.L_d-emag.L_q).*-x(1).^2.*0.5.*sin(2.*x(2)) + emag.Psi_PM.*-x(1).*sin(x(2)));
    % Starting value
    x0 = [0 0];
    % Bound Constraints
    lb = [0 -pi];
    ub = [opt.i_1max 0];
else
    % Minimization function 
    fun = @(x) -1.5.*rated.p.*((emag.L_d-emag.L_q).*x(1).^2.*0.5.*sin(2.*x(2)) + emag.Psi_PM.*x(1).*sin(x(2)));
    % Starting value
    x0 = [0 0];
    % Bound Constraints
    lb = [0 0];
    ub = [opt.i_1max pi];
end
% Linear Inequality Constraint
A = [];
b = [];
% Linear Equality Constraint
Aeq = [];
beq = [];
% Nonlinear Constraints
nonlcon = @nonlcon_Optimierung_M;
% Option for Optimization
options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point','MaxFunctionEvaluations',3e4,'MaxIterations',2000,...
                            'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-6);

% Memory allokation
i_1 = zeros(length(omega_vec(1,:)),1);
beta = zeros(length(omega_vec(1,:)),1);

tic
for j = 1:length(omega_vec(1,:))
    omega_k = omega_vec(1,j);
    % Minimization with start value [0,0]
    [x_sol,~,exitflag,~] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, ...
    @(x) nonlcon(x, opt.u_1max, emag.L_d, emag.L_q, emag.R_1, omega_k, emag.Psi_PM), options);
    if(exitflag~=0 && exitflag~=1 && exitflag~=2)
        % Current is set inf if no result is found
        i_1(j) = inf;
        beta(j) = inf;
    else
        % Assign the values when result is found
        i_1(j) = x_sol(1);
        beta(j) = x_sol(2);
    end
end
toc

end

function [c,ceq] = nonlcon_Optimierung_M(x, u_1max, L_d, L_q, R_1, omega_k, Psi_PM)
%nonlcon_fun Definition of nonlinear constraints
%   Nonlinear Inequality Constraint: voltage
%   Nonlinear Equality Constraint: none

    c = (R_1*x(1)*cos(x(2)) - omega_k*L_q*x(1)*sin(x(2)))^2 + (R_1*x(1)*sin(x(2)) + omega_k*L_d*x(1)*cos(x(2)) + omega_k*Psi_PM)^2 - u_1max^2;
    ceq = [];
end

%% F.2) Optimierung_i - Current optimization
% #########################################################################
% #   F.2) OPTIMIERUNG_I - CURRENT OPTIMIZATION                           #
% #########################################################################

function [i_1dsol, i_1qsol] = Optimierung_i(rated, emag, opt, omega_mesh, M_mesh, Generator)
% Optimization Minimizes the d and q components of the current under
% constraints and limits
% Notation: i_1d = x(1), i_1q = x(2)
% Bound (motor): -i_1max<=i_1d<=0 (for L_d<L_q) or 0<=i1_d<=i_1max (for L_d>L_q), 0<=i_1q<=i_1max
% Linear Inequality Constraint: none
% Linear Equality Constraint: none
% Nonlinear Constraints: see function nonlcon_Optimierung_i
% For reasons of calculation time, from column 2 or line 2 values of
% of previous currents are used as start values.

% Minimization function 
fun = @(x) x(1)^2 + x(2)^2;

if(Generator)
    % Start value
    x0 = [0 0];
    % Bound Constraints
    if(emag.L_d<=emag.L_q)
        lb = [-opt.i_1max -opt.i_1max];
        ub = [0 0];
    else
        lb = [0 -opt.i_1max];
        ub = [opt.i_1max 0];
    end
else
    % Start value
    x0 = [0 0];
    % Bound Constraints
    if(emag.L_d<=emag.L_q)
        lb = [-opt.i_1max 0];
        ub = [0 opt.i_1max];
    else
        lb = [0 0];
        ub = [opt.i_1max opt.i_1max];
    end
end
% Linear Inequality Constraint
A = [];
b = [];
% Linear Equality Constraint
Aeq = [];
beq = [];
% Nonlinear Constraints
nonlcon = @nonlcon_Optimierung_i;
% Options for Optimization
options = optimoptions(@fmincon,'Display','off','Algorithm','sqp','MaxFunctionEvaluations',3e4,'MaxIterations',400,...
                            'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-6);

% Memory allokation
i_1dsol = zeros(size(M_mesh));
i_1qsol = zeros(size(M_mesh));

tic
for j = 1:length(omega_mesh(1,:))
    omega_k = omega_mesh(1,j);
    for k = 1:length(M_mesh(:,1))
        M = M_mesh(k,j);
        if(isnan(M))
            % Sorting out the points that cannot be approached (full load characteristic)
            i_1dsol(k,j) = NaN;
            i_1qsol(k,j) = NaN;
        else
            if(j>1 && ~isnan(i_1dsol(k,j-1)) && ~isnan(i_1qsol(k,j-1)))
                if(k>1 && ~isnan(i_1dsol(k-1,j)))
                    % Minimization with initial values from previous calculation
                    [x_sol,~,exitflag,~] = fmincon(fun, [i_1dsol(k-1,j) i_1qsol(k,j-1)], A, b, Aeq, beq, lb, ub, ...
                    @(x) nonlcon(x, M, opt.u_1max, opt.i_1max, emag.L_d, emag.L_q, emag.R_1, omega_k, emag.Psi_PM, rated.p), options);
                else
                    % Minimization with initial values from previous calculation
                    [x_sol,~,exitflag,~] = fmincon(fun, [i_1dsol(k,j-1) i_1qsol(k,j-1)], A, b, Aeq, beq, lb, ub, ...
                        @(x) nonlcon(x, M, opt.u_1max, opt.i_1max, emag.L_d, emag.L_q, emag.R_1, omega_k, emag.Psi_PM, rated.p), options);
                end
            else
                % Minimization with start value [0,0].
                [x_sol,~,exitflag,~] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, ...
                @(x) nonlcon(x, M, opt.u_1max, opt.i_1max, emag.L_d, emag.L_q, emag.R_1, omega_k, emag.Psi_PM, rated.p), options);
            end
            % values < tol set to 0
            if(abs(x_sol(1))<1e-2)
                x_sol(1)=0;
            end
            if(abs(x_sol(2))<1e-2)
                x_sol(2)=0;
            end
            if(exitflag~=0 && exitflag~=1 && exitflag~=2)
                % Minimization with startvalue [0,0]
                [x_sol,~,exitflag,~] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, ...
                @(x) nonlcon(x, M, opt.u_1max, opt.i_1max, emag.L_d, emag.L_q, emag.R_1, omega_k, emag.Psi_PM, rated.p), options);
                if(exitflag~=0 && exitflag~=1 && exitflag~=2)
                    % Current is set inf if no result is found
                    i_1dsol(k,j) = inf;
                    i_1qsol(k,j) = inf;
                     warning(['fmincon exitflag:', num2str(exitflag), ', check function Optimierung_i'])
                else
                    % Current is set inf if no result is found
                    i_1dsol(k,j) = x_sol(1);
                    i_1qsol(k,j) = x_sol(2);
                end
            else
                % Assign the values when result is found
                i_1dsol(k,j) = x_sol(1);
                i_1qsol(k,j) = x_sol(2);
            end
        end
    end
end
toc

end

function [c,ceq] = nonlcon_Optimierung_i(x, M_ref, u_1max, i_1max, L_d, L_q, R_1, omega_k, Psi_PM, p)
%nonlcon_fun Definition of the nonlinear constraints
% Nonlinear Inequality Constraint: current limit, voltage limit
% Nonlinear Equality Constraint: Moment Reference
    
    c = [x(1)^2 + x(2)^2 - i_1max^2;...
        (R_1*x(1)-omega_k*L_q*x(2))^2 + (R_1*x(2)+omega_k*L_d*x(1)+omega_k*Psi_PM)^2 - u_1max^2];
    ceq = 1.5*p*(L_d-L_q)*x(2)*x(1) + 1.5*p*Psi_PM*x(2) - M_ref;
end