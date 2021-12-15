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
% B.1) Engine model
% B.2) Generator model
% C) Loss calculation
% D) Map generation
% E) Postprocessing
% F) Auxiliary functions
% F.1) Optimization_M
% F.2) Optimization_i

function [Analyse] = Analysis_ASM(handles)

%% A) Preprocessing
% #########################################################################
% #   A) PREPROCESSING                                                    #
% #########################################################################

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


betr.omega_k_max = (opt.n_max/60)*(2*pi)*rated.p;


betr.omega_k_vec = linspace(0, betr.omega_k_max, opt.n_tics);

%% B.1) Motormodel
% #########################################################################
% #   B.1) MOTORMODELL                                                    #
% #########################################################################

[ctrl.mot.i_1, ctrl.mot.beta] = Optimierung_M(rated, emag, opt, betr.omega_k_vec, 0);
betr.mot.M_max_vec = 1.5.*rated.p.*((emag.L_12.^2./emag.L_22).*ctrl.mot.i_1.^2.*0.5.*sin(2.*ctrl.mot.beta));
betr.mot.M_max_vec = betr.mot.M_max_vec';

tol = 1e-3;

betr.mot.M_vec = linspace(max(betr.mot.M_max_vec)-tol, 0, opt.M_tics);
clear tol

[betr.mot.omega_k_mesh, betr.mot.M_max_mesh] = meshgrid(betr.omega_k_vec, betr.mot.M_vec);

for i = 1:length(betr.mot.omega_k_mesh(1,:))
    betr.mot.M_max_mesh(betr.mot.M_max_mesh(:,i)>betr.mot.M_max_vec(i),i) = NaN;
end


[ctrl.mot.i_1d_mesh, ctrl.mot.i_1q_mesh] = Optimierung_i(rated, emag, opt, betr.mot.omega_k_mesh, betr.mot.M_max_mesh, 0);


ctrl.mot.i_2d_mesh = zeros(size(betr.mot.M_max_mesh));
ctrl.mot.i_2q_mesh = -(emag.L_12./emag.L_22) .* ctrl.mot.i_1q_mesh;


ctrl.mot.u_1d_mesh = emag.R_1.*ctrl.mot.i_1d_mesh - betr.mot.omega_k_mesh.*emag.L_11.*emag.sigma.*ctrl.mot.i_1q_mesh;
ctrl.mot.u_1q_mesh = emag.R_1.*ctrl.mot.i_1q_mesh + betr.mot.omega_k_mesh.*emag.L_11.*ctrl.mot.i_1d_mesh;

ctrl.mot.P_el_mesh = 1.5.*(ctrl.mot.i_1d_mesh.*ctrl.mot.u_1d_mesh + ctrl.mot.i_1q_mesh.*ctrl.mot.u_1q_mesh);


ctrl.mot.i_1d_mesh(end,1) = 0;
ctrl.mot.i_1q_mesh(end,1) = 0;
var1 = ctrl.mot.i_1d_mesh(end,:);
var2 = ctrl.mot.i_1q_mesh(end,:);
if(all(ctrl.mot.i_1d_mesh(end,:)<1e-10) && all(ctrl.mot.i_1q_mesh(end,:)<1e-10)) 
    ctrl.mot.i_1d_mesh(end,:) = 1e-10;
    ctrl.mot.i_1q_mesh(end,:) = 1e-10;
end
betr.mot.omega_el_mesh = betr.mot.omega_k_mesh - ((emag.R_2./emag.L_22).*(ctrl.mot.i_1q_mesh./ctrl.mot.i_1d_mesh)); 
betr.mot.omega_el_mesh(end,1) = 0;
ctrl.mot.i_1d_mesh(end,:) = var1;
ctrl.mot.i_1q_mesh(end,:) = var2;
betr.mot.f_el_mesh = betr.mot.omega_el_mesh ./ (2.*pi);  
betr.mot.omega_m_mesh = (betr.mot.omega_el_mesh ./ rated.p); 
betr.mot.n_m_mesh = (betr.mot.omega_m_mesh ./ (2.*pi)) .* 60; 
for i = 1:length(betr.mot.n_m_mesh(1,:))
    var = betr.mot.n_m_mesh(~isnan(betr.mot.n_m_mesh(:,i)),i);
    betr.mot.n_m_vec(i) = var(1);
end
clear var1 var2 var

%% B.2) Generatormodel
% #########################################################################
% #   B.2) GENERATORMODEL                                                 #
% #########################################################################

if(opt.Generator)
    
    [ctrl.gen.i_1, ctrl.gen.beta] = Optimierung_M(rated, emag, opt, betr.omega_k_vec, 1);
    betr.gen.M_max_vec = 1.5.*rated.p.*((emag.L_12.^2./emag.L_22).*ctrl.gen.i_1.^2.*0.5.*sin(2.*ctrl.gen.beta));
    betr.gen.M_max_vec = betr.gen.M_max_vec';
    
   
    tol = 1e-3;
    
    
    betr.gen.M_vec = linspace(0, min(betr.gen.M_max_vec)+tol, opt.M_tics);
    clear tol

    
    [betr.gen.omega_k_mesh, betr.gen.M_max_mesh] = meshgrid(betr.omega_k_vec, betr.gen.M_vec);
    
    for i = 1:length(betr.gen.omega_k_mesh(1,:))
        betr.gen.M_max_mesh(betr.gen.M_max_mesh(:,i)<betr.gen.M_max_vec(i),i) = NaN;
    end
    
    
    [ctrl.gen.i_1d_mesh, ctrl.gen.i_1q_mesh] = Optimierung_i(rated, emag, opt, betr.gen.omega_k_mesh, betr.gen.M_max_mesh, 1);
    
    
    ctrl.gen.i_2d_mesh = zeros(size(betr.gen.M_max_mesh));
    ctrl.gen.i_2q_mesh = -(emag.L_12./emag.L_22) .* ctrl.gen.i_1q_mesh;


    ctrl.gen.u_1d_mesh = emag.R_1.*ctrl.gen.i_1d_mesh - betr.gen.omega_k_mesh.*emag.L_11.*emag.sigma.*ctrl.gen.i_1q_mesh;
    ctrl.gen.u_1q_mesh = emag.R_1.*ctrl.gen.i_1q_mesh + betr.gen.omega_k_mesh.*emag.L_11.*ctrl.gen.i_1d_mesh;

    
    ctrl.gen.P_el_mesh = 1.5.*(ctrl.gen.i_1d_mesh.*ctrl.gen.u_1d_mesh + ctrl.gen.i_1q_mesh.*ctrl.gen.u_1q_mesh);

    
    ctrl.gen.i_1d_mesh(end,1) = 0;
    ctrl.gen.i_1q_mesh(end,1) = 0;
    var1 = ctrl.gen.i_1d_mesh(1,:);
    var2 = ctrl.gen.i_1q_mesh(1,:);
    if(all(ctrl.gen.i_1d_mesh(1,:)<1e-10) && all(ctrl.gen.i_1q_mesh(1,:)<1e-10)) 
        ctrl.gen.i_1d_mesh(1,:) = 1e-10;
        ctrl.gen.i_1q_mesh(1,:) = 1e-10;
    end
    betr.gen.omega_el_mesh = betr.gen.omega_k_mesh - ((emag.R_2./emag.L_22).*(ctrl.gen.i_1q_mesh./ctrl.gen.i_1d_mesh)); 
    betr.gen.omega_el_mesh(1,1) = 0;
    ctrl.gen.i_1d_mesh(1,:) = var1;
    ctrl.gen.i_1q_mesh(1,:) = var2;
    betr.gen.f_el_mesh = betr.gen.omega_el_mesh ./ (2.*pi);  
    betr.gen.omega_m_mesh = (betr.gen.omega_el_mesh ./ rated.p); 
    betr.gen.n_m_mesh = (betr.gen.omega_m_mesh ./ (2.*pi)) .* 60; 
    for i = 1:length(betr.gen.n_m_mesh(1,:))
        var = betr.gen.n_m_mesh(~isnan(betr.gen.n_m_mesh(:,i)),i);
        betr.gen.n_m_vec(i) = var(end);
    end
    clear var1 var2 var
    
    betr.M_vec = [betr.mot.M_vec betr.gen.M_vec(:,2:end)];
    betr.omega_k_mesh = [betr.mot.omega_k_mesh; betr.gen.omega_k_mesh(2:end,:)];
    betr.M_max_mesh = [betr.mot.M_max_mesh; betr.gen.M_max_mesh(2:end,:)];
    ctrl.i_1d_mesh = [ctrl.mot.i_1d_mesh; ctrl.gen.i_1d_mesh(2:end,:)];
    ctrl.i_1q_mesh = [ctrl.mot.i_1q_mesh; ctrl.gen.i_1q_mesh(2:end,:)];
    ctrl.i_2d_mesh = [ctrl.mot.i_2d_mesh; ctrl.gen.i_2d_mesh(2:end,:)];
    ctrl.i_2q_mesh = [ctrl.mot.i_2q_mesh; ctrl.gen.i_2q_mesh(2:end,:)];
    ctrl.u_1d_mesh = [ctrl.mot.u_1d_mesh; ctrl.gen.u_1d_mesh(2:end,:)];
    ctrl.u_1q_mesh = [ctrl.mot.u_1q_mesh; ctrl.gen.u_1q_mesh(2:end,:)];
    ctrl.P_el_mesh = [ctrl.mot.P_el_mesh; ctrl.gen.P_el_mesh(2:end,:)];
    betr.omega_el_mesh = [betr.mot.omega_el_mesh; betr.gen.omega_el_mesh(2:end,:)];
    betr.f_el_mesh = [betr.mot.f_el_mesh; betr.gen.f_el_mesh(2:end,:)];
    betr.omega_m_mesh = [betr.mot.omega_m_mesh; betr.gen.omega_m_mesh(2:end,:)];
    betr.n_m_mesh = [betr.mot.n_m_mesh; betr.gen.n_m_mesh(2:end,:)];
    % betr.n_m_vec = [betr.mot.n_m_vec betr.gen.n_m_vec(:,2:end)];
else
   
    % ctrl.i_1 = ctrl.mot.i_1;
    % ctrl.beta = ctrl.mot.beta;
    % betr.M_max_vec = betr.mot.M_max_vec;
    betr.M_vec = betr.mot.M_vec;
    betr.omega_k_mesh = betr.mot.omega_k_mesh;
    betr.M_max_mesh = betr.mot.M_max_mesh;
    ctrl.i_1d_mesh = ctrl.mot.i_1d_mesh;
    ctrl.i_1q_mesh = ctrl.mot.i_1q_mesh;
    ctrl.i_2d_mesh = ctrl.mot.i_2d_mesh;
    ctrl.i_2q_mesh = ctrl.mot.i_2q_mesh;
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


% (i_d_mesh.^2 + i_q_mesh.^2) sind peak Werte
if(opt.P_vw==1)
    loss.vw.P_1vw_mesh = 1.5 * emag.R_1 * (ctrl.i_1d_mesh.^2 + ctrl.i_1q_mesh.^2); %Stator
    loss.vw.P_2vw_mesh = 1.5 * emag.R_2 * (ctrl.i_2d_mesh.^2 + ctrl.i_2q_mesh.^2); %Rotor
    loss.P_vw_mesh = loss.vw.P_1vw_mesh + loss.vw.P_2vw_mesh;
else
    loss.P_vw_mesh = zeros(size(ctrl.i_1d_mesh));
end

% P_vu [W]
if(opt.P_vu==1)
    if(strcmp(opt.P_vu_Modell,'Modellansatz Jordan'))
       
        loss.vu.w_1Fe = opt_Entwurf.Stator_Eisenmaterial.p_vFe(opt_Entwurf.Stator_Eisenmaterial.p_vFe_B_vec==1,:) ./ opt_Entwurf.Stator_Eisenmaterial.p_vFe_f_vec;
        [loss.vu.fitresult_1, ~] = fit(opt_Entwurf.Stator_Eisenmaterial.p_vFe_f_vec', loss.vu.w_1Fe', 'poly1');
        loss.vu.fitresult_data_1 = coeffvalues(loss.vu.fitresult_1);
        loss.vu.sigma_1h_s = loss.vu.fitresult_data_1(2);
        loss.vu.sigma_1w_s = loss.vu.fitresult_data_1(1);
        loss.vu.sigma_1h = loss.vu.sigma_1h_s * 50;
        loss.vu.sigma_1w = loss.vu.sigma_1w_s * 50^2;
        
       
        loss.vu.w_2Fe = opt_Entwurf.Stator_Eisenmaterial.p_vFe(opt_Entwurf.Rotor_Eisenmaterial.p_vFe_B_vec==1,:) ./ opt_Entwurf.Stator_Eisenmaterial.p_vFe_f_vec;
        [loss.vu.fitresult_2, ~] = fit(opt_Entwurf.Stator_Eisenmaterial.p_vFe_f_vec', loss.vu.w_2Fe', 'poly1');
        loss.vu.fitresult_data_2 = coeffvalues(loss.vu.fitresult_2);
        loss.vu.sigma_2h_s = loss.vu.fitresult_data_2(2);
        loss.vu.sigma_2w_s = loss.vu.fitresult_data_2(1);
        loss.vu.sigma_2h = loss.vu.sigma_2h_s * 50;
        loss.vu.sigma_2w = loss.vu.sigma_2w_s * 50^2;
        
       
        loss.vu.c_h = 1;
        loss.vu.c_w = 1;
        
        
        B_1r = emag.B_1r * ones(size(betr.f_el_mesh));
        B_1z_m = emag.B_1z_m * ones(size(betr.f_el_mesh));
        loss.vu.p_1r = (loss.vu.sigma_1h.*loss.vu.c_h.*(betr.f_el_mesh./50) + loss.vu.sigma_1w.*loss.vu.c_w.*(betr.f_el_mesh./50).^2) .* (B_1r./1).^2;
        loss.vu.p_1z = (loss.vu.sigma_1h.*loss.vu.c_h.*(betr.f_el_mesh./50) + loss.vu.sigma_1w.*loss.vu.c_w.*(betr.f_el_mesh./50).^2) .* (B_1z_m./1).^2;
        
        
        loss.vu.k_u_r = 1.7; 
        loss.vu.k_u_z = 1.8; 

        
        loss.vu.P_1r = loss.vu.k_u_r .* loss.vu.p_1r .* opt_Entwurf.Stator_Eisenmaterial.rho_Fe .* geo.Vo_1r;
        loss.vu.P_1z = loss.vu.k_u_z .* loss.vu.p_1z .* opt_Entwurf.Stator_Eisenmaterial.rho_Fe .* geo.Vo_1z;
        
       
        loss.P_vu_mesh = loss.vu.P_1r + loss.vu.P_1z;
    else
        error('Ungueltige Eingabe bei Variable "opt.P_vu_Modell"');
    end
else
    loss.P_vu_mesh = zeros(size(ctrl.i_1d_mesh));
end


if(opt.P_vme==1)
    % k_rb [Ws^2/m^4]
    loss.vme.k_rb = 10;

    % [m/s]
    loss.vme.v_2 = (geo.D_2a ./ 2) .* betr.omega_m_mesh;

    loss.P_vme_mesh = loss.vme.k_rb .* geo.D_2a .* (geo.l_i + 0.8^3 .* 0.6 .* geo.tau_1p) .* loss.vme.v_2.^2;
else
	loss.P_vme_mesh = zeros(size(ctrl.i_1d_mesh));
end

%  P_zus [W]

if(opt.P_vzus==1)
    loss.P_vzus_mesh = abs(ctrl.P_el_mesh) .* (0.025 - 0.005*log(rated.P_N*1e-3));
else
    loss.P_vzus_mesh = zeros(size(ctrl.i_1d_mesh));
end

%P_ges [W]
loss.P_vges_mesh = loss.P_vw_mesh + loss.P_vu_mesh + loss.P_vme_mesh + loss.P_vzus_mesh;

%% D) Map creation
% #########################################################################
% #   D) MAP CREATION                                                     #
% #########################################################################


betr.mot.P_mech_mesh = betr.mot.M_max_mesh .* betr.mot.omega_m_mesh;


[var1, ~] = size(betr.mot.P_mech_mesh);
eta.mot.eta_vw_mesh_alt = betr.mot.P_mech_mesh ./ ctrl.mot.P_el_mesh;
eta.mot.eta_vw_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vw_mesh(1:var1,1:end));
eta.mot.eta_fe_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vu_mesh(1:var1,1:end));
eta.mot.eta_vme_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vme_mesh(1:var1,1:end));
eta.mot.eta_zus_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vzus_mesh(1:var1,1:end));
eta.mot.eta_ges_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vges_mesh(1:var1,1:end));

if(opt.Generator)
   
    betr.gen.P_mech_mesh = betr.gen.M_max_mesh .* betr.gen.omega_m_mesh;
    
    
    eta.gen.eta_vw_mesh_alt = ctrl.gen.P_el_mesh ./ betr.gen.P_mech_mesh;
    eta.gen.eta_vw_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vw_mesh(var1:end,1:end));
    eta.gen.eta_fe_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vu_mesh(var1:end,1:end));
    eta.gen.eta_vme_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vme_mesh(var1:end,1:end));
    eta.gen.eta_zus_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vzus_mesh(var1:end,1:end));
    eta.gen.eta_ges_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vges_mesh(var1:end,1:end));
    
  
    eta.eta_vw_mesh = [eta.mot.eta_vw_mesh; eta.gen.eta_vw_mesh(2:end,:)];
    eta.eta_fe_mesh = [eta.mot.eta_fe_mesh; eta.gen.eta_fe_mesh(2:end,:)];
    eta.eta_vme_mesh = [eta.mot.eta_vme_mesh; eta.gen.eta_vme_mesh(2:end,:)];
    eta.eta_zus_mesh = [eta.mot.eta_zus_mesh; eta.gen.eta_zus_mesh(2:end,:)];
    eta.eta_ges_mesh = [eta.mot.eta_ges_mesh; eta.gen.eta_ges_mesh(2:end,:)];
else
    
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


% betr.n_m_max = max([max(max(betr.mot.n_m_mesh)); max(max(betr.gen.n_m_mesh))]);
% betr.n_m_min = max([min(betr.mot.n_m_mesh(1:end-1,1)); min(betr.gen.n_m_mesh(2:end,1))]);
% betr.n_m_vec_spaced = linspace(betr.n_m_min, betr.n_m_max, opt.n_tics);
betr.n_m_vec_spaced = linspace(0, opt.n_max, opt.n_tics);
[betr.n_m_mesh_spaced, ~] = meshgrid(betr.n_m_vec_spaced, betr.M_vec);
for i = 1:length(betr.n_m_mesh(:,1))
    var = interp1(betr.n_m_mesh(i,~isnan(ctrl.i_1d_mesh(i,:))), ctrl.i_1d_mesh(i,~isnan(ctrl.i_1d_mesh(i,:))), betr.n_m_mesh_spaced(i,~isnan(ctrl.i_1d_mesh(i,:))), 'linear','extrap');
    ctrl.i_d_mesh_spaced(i,:) = [var NaN(1,length(betr.n_m_mesh(i,isnan(ctrl.i_1d_mesh(i,:)))))];
    var = interp1(betr.n_m_mesh(i,~isnan(ctrl.i_1q_mesh(i,:))), ctrl.i_1q_mesh(i,~isnan(ctrl.i_1q_mesh(i,:))), betr.n_m_mesh_spaced(i,~isnan(ctrl.i_1q_mesh(i,:))), 'linear','extrap');
    ctrl.i_q_mesh_spaced(i,:) = [var NaN(1,length(betr.n_m_mesh(i,isnan(ctrl.i_1q_mesh(i,:)))))];
    var = interp1(betr.n_m_mesh(i,~isnan(ctrl.u_1d_mesh(i,:))), ctrl.u_1d_mesh(i,~isnan(ctrl.u_1d_mesh(i,:))), betr.n_m_mesh_spaced(i,~isnan(ctrl.u_1d_mesh(i,:))), 'linear','extrap');
    ctrl.u_d_mesh_spaced(i,:) = [var NaN(1,length(betr.n_m_mesh(i,isnan(ctrl.u_1d_mesh(i,:)))))];
    var = interp1(betr.n_m_mesh(i,~isnan(ctrl.u_1q_mesh(i,:))), ctrl.u_1q_mesh(i,~isnan(ctrl.u_1q_mesh(i,:))), betr.n_m_mesh_spaced(i,~isnan(ctrl.u_1q_mesh(i,:))), 'linear','extrap');
    ctrl.u_q_mesh_spaced(i,:) = [var NaN(1,length(betr.n_m_mesh(i,isnan(ctrl.u_1q_mesh(i,:)))))];
end
clear i var

opt.Locked = 1;


ctrl = orderfields(ctrl);
ctrl.mot = orderfields(ctrl.mot);
if(opt.Generator)
    ctrl.gen = orderfields(ctrl.gen);
    betr.gen = orderfields(betr.gen);
end
loss = orderfields(loss);
if(isfield(loss,'vw'))
    loss.vw = orderfields(loss.vw);
end
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


Analyse.Momentensteuerung = ctrl;
Analyse.Verluste = loss;
Analyse.Betriebsdaten = betr;
Analyse.Wirkungsgrad = eta;
Analyse.Optionen = opt;
Analyse.EMAG = emag;
Analyse.Geometrie = geo;
Analyse.Optionen = opt;


%save(['3_Results/',handles.Entwurf.Optionen.folder_id ,'/2_Analyse/Analyse_',handles.Entwurf.Optionen.file_id  ,'.mat'],'Analyse');

end

%% F) Help functions
% #########################################################################
% #   F) HELPFUNCTIONS                                                  #
% #########################################################################

%% F.1) Optimization_M
% #########################################################################
% #   F.1) OPTIMIZATION_M                                                 #
% #########################################################################

function [i_1, beta] = Optimierung_M(rated, emag, opt, omega_vec, Generator)


if(Generator)
    fun = @(x) -1.5.*rated.p.*((emag.L_12.^2./emag.L_22).*-x(1).^2.*0.5.*sin(2.*x(2)));
    
    x0 = [0 0];
    lb = [0 -pi];
    ub = [opt.i_1max 0];
else
   
    fun = @(x) -1.5.*rated.p.*((emag.L_12.^2./emag.L_22).*x(1).^2.*0.5.*sin(2.*x(2)));
    
   
    x0 = [0 0];
   
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

options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point','MaxFunctionEvaluations',3e4,'MaxIterations',2000,...
                            'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-6);


i_1 = zeros(length(omega_vec(1,:)),1);
beta = zeros(length(omega_vec(1,:)),1);

tic
for j = 1:length(omega_vec(1,:))
    omega_k = omega_vec(1,j);
  
    [x_sol,~,exitflag,~] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, ...
    @(x) nonlcon(x, opt.u_1max, emag.L_11, emag.sigma, emag.R_1, omega_k), options);
    if(exitflag~=0 && exitflag~=1 && exitflag~=2)
        
        i_1(j) = inf;
        beta(j) = inf;
    else
       
        i_1(j) = x_sol(1);
        beta(j) = x_sol(2);
    end
end
toc

end

function [c,ceq] = nonlcon_Optimierung_M(x, u_1max, L_11, sigma, R_1, omega_k)
%nonlcon_fun Definition der nichtlinearen Nebenbedingungen
%   Nonlinear Inequality Constraint: Spannungsgrenze
%   Nonlinear Equality Constraint: keine
    
    c = (R_1*x(1)*cos(x(2)) - omega_k*L_11*sigma*x(1)*sin(x(2)))^2 + (R_1*x(1)*sin(x(2)) + omega_k*L_11*x(1)*cos(x(2)))^2 - u_1max^2;
    ceq = [];
end

%% F.2) Optimization_i
% #########################################################################
% #   F.2) OPTIMIZATION_I                                                 #
% #########################################################################

function [i_1dsol, i_1qsol] = Optimierung_i(rated, emag, opt, omega_mesh, M_mesh, Generator)

fun = @(x) x(1)^2 + x(2)^2;

if(Generator)
 
    x0 = [opt.i_1max -opt.i_1max];

    lb = [0 -opt.i_1max];
    ub = [opt.i_1max 0];
else
    
    x0 = [opt.i_1max opt.i_1max];
  
    lb = [0 0];
    ub = [opt.i_1max opt.i_1max];
end

% Linear Inequality Constraint
A = [];
b = [];
% Linear Equality Constraint
Aeq = [];
beq = [];
% Nonlinear Constraints
nonlcon = @nonlcon_Optimierung_i;

options = optimoptions(@fmincon,'Display','off','Algorithm','sqp','MaxFunctionEvaluations',3e4,'MaxIterations',400,...
                            'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-6);


i_1dsol = zeros(size(M_mesh));
i_1qsol = zeros(size(M_mesh));

tic
for j = 1:length(omega_mesh(1,:))
    omega_k = omega_mesh(1,j);
    for k = 1:length(M_mesh(:,1))
        M = M_mesh(k,j);
        if(isnan(M))
            
            i_1dsol(k,j) = NaN;
            i_1qsol(k,j) = NaN;
        else
            if(j>1 && ~isnan(i_1dsol(k,j-1)) && ~isnan(i_1qsol(k,j-1)))
                if(k>1 && ~isnan(i_1dsol(k-1,j)))
                    
                    [x_sol,~,exitflag,~] = fmincon(fun, [i_1dsol(k-1,j) i_1qsol(k,j-1)], A, b, Aeq, beq, lb, ub, ...
                    @(x) nonlcon(x, M, opt.u_1max, opt.i_1max, emag.L_11, emag.L_12, emag.L_22, emag.sigma, emag.R_1, omega_k, rated.p), options);
                else
                    
                    [x_sol,~,exitflag,~] = fmincon(fun, [i_1dsol(k,j-1) i_1qsol(k,j-1)], A, b, Aeq, beq, lb, ub, ...
                        @(x) nonlcon(x, M, opt.u_1max, opt.i_1max, emag.L_11, emag.L_12, emag.L_22, emag.sigma, emag.R_1, omega_k, rated.p), options);
                end
            else
                
                [x_sol,~,exitflag,~] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, ...
                @(x) nonlcon(x, M, opt.u_1max, opt.i_1max, emag.L_11, emag.L_12, emag.L_22, emag.sigma, emag.R_1, omega_k, rated.p), options);
            end
            )
            if(abs(x_sol(1))<1e-2)
                x_sol(1)=0;
            end
            if(abs(x_sol(2))<1e-2)
                x_sol(2)=0;
            end
            if(exitflag~=0 && exitflag~=1 && exitflag~=2)
                
                [x_sol,~,exitflag,~] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, ...
                @(x) nonlcon(x, M, opt.u_1max, opt.i_1max, emag.L_11, emag.L_12, emag.L_22, emag.sigma, emag.R_1, omega_k, rated.p), options);
                if(exitflag~=0 && exitflag~=1 && exitflag~=2)
                    
                    i_1dsol(k,j) = inf;
                    i_1qsol(k,j) = inf;
                    warning(['fmincon exitflag:', num2str(exitflag), ', check function Optimierung_i'])
                else
                    
                    i_1dsol(k,j) = x_sol(1);
                    i_1qsol(k,j) = x_sol(2);
                end
            else
                
                i_1dsol(k,j) = x_sol(1);
                i_1qsol(k,j) = x_sol(2);
            end
        end
    end
end
toc

end

function [c,ceq] = nonlcon_Optimierung_i(x, M_ref, u_1max, i_1max, L_11, L_12, L_22, sigma, R_1, omega_k, p)
%nonlcon_fun Definition der nichtlinearen Nebenbedingungen
%   Nonlinear Inequality Constraint: Stromgrenze, Spannungsgrenze
%   Nonlinear Equality Constraint: Momentenreferenz
    
    c = [x(1)^2 + x(2)^2 - i_1max^2;...
        ((R_1*x(1)) - (omega_k*L_11*sigma*x(2)))^2 + ((R_1*x(2)) + (omega_k*L_11*x(1)))^2 - u_1max^2];
    ceq = 1.5*p*(L_12^2/L_22)*x(2)*x(1) - M_ref;
end