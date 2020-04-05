% -------------------------------------------------------------------------
% TU Muenchen - Lehrstuhl fuer Fahrzeugtechnik (FTM)
% -------------------------------------------------------------------------
% Modell fuer den Entwurf und die Analyse einer PMSM oder ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor: Svenja Kalt (kalt@ftm.mw.tum.de)
%        Jonathan Erhard
% -------------------------------------------------------------------------

%% Inhaltsverzeichnis
% A) Preprocessing
% B) Momentensteuerung
%   B.1) Motormodell
%   B.2) Generatormodell
% C) Verlustberechnung
% D) Kennfelderstellung
% E) Postprocessing
% F) Hilfsfunktionen
%   F.1) Optimierung_M
%   F.2) Optimierung_i

function [Analyse] = Analyse_ASM(handles)
%Analyse Diese Funktion fuehrt die Analyse einer ASM durch.
%   In der Funktion Analyse_ASM werden Kennfelder fuer die vorher
%   ausgelegte Maschine berechnet. Dazu werden zuerst die Stroeme und
%   Spannungen mittels eines Motormodells berechnet. Dann werden die
%   einzelnen Verlustanteile ermittelt. Zuletzt entsteht das
%   Wirkungsgradkennfeld. Hinweis: Es erfolgt keine thermische
%   Ueberpruefung o.ae. Der Nutzer ist deswegen angehalten die Maschine nur
%   mit realistischen Angaben zu speisen (i_max und u_max).
%   V/A: lineare & stationaere Betrachtung:
%   - keine Saettigung
%   - keine Kreuzkopplung der Induktivitaeten
%   - Induktivitaeten sind nicht stromabhaengig
%   - keine Dynamik

%% A) Preprocessing
% #########################################################################
% #   A) PREPROCESSING                                                    #
% #########################################################################

% Umspeichern der Variablen zur einfacheren Nutzung
const = handles.Entwurf.Konstanten;
emag = handles.Entwurf.EMAG;
geo = handles.Entwurf.Geometrie;
wick = handles.Entwurf.Wicklung;

rated = handles.Entwurf.Bemessungswerte;
richt = handles.Entwurf.Richtwerte;
opt_Entwurf = handles.Entwurf.Optionen;

opt = handles.opt;

%% B) Momentensteuerung
% #########################################################################
% #   B) MOMENTENSTEUERUNG                                                #
% #########################################################################

% Drehzahlvektor
% V/A: omega_k_max ungefaehr gleich omega_m_max fuer max. Drehzahl im Kennfeld
betr.omega_k_max = (opt.n_max/60)*(2*pi)*rated.p;

% Erstellen Drehzahlvektor
betr.omega_k_vec = linspace(0, betr.omega_k_max, opt.n_tics); % elektrische Winkelgeschwindigkeit (k-Koordinatensystem)

%% B.1) Motormodell
% #########################################################################
% #   B.1) MOTORMODELL                                                    #
% #########################################################################

% Berechnung der Volllastkennlinie (maximales Moment fuer Drehzahlvektor berechnen)
% Maximierung des Moments unter Nebenbedingungen (Betriebsbereich Motor)
[ctrl.mot.i_1, ctrl.mot.beta] = Optimierung_M(rated, emag, opt, betr.omega_k_vec, 0);
betr.mot.M_max_vec = 1.5.*rated.p.*((emag.L_12.^2./emag.L_22).*ctrl.mot.i_1.^2.*0.5.*sin(2.*ctrl.mot.beta));
betr.mot.M_max_vec = betr.mot.M_max_vec';

% Toleranz für max Moment (Grund: numerische Ungenauigkeit)
tol = 1e-3;

% Erstellen Momentenvektor
betr.mot.M_vec = linspace(max(betr.mot.M_max_vec)-tol, 0, opt.M_tics);
clear tol

% Erstellen Drehzahl-/Drehmomentgitter
[betr.mot.omega_k_mesh, betr.mot.M_max_mesh] = meshgrid(betr.omega_k_vec, betr.mot.M_vec);

% Anpassung Kennfeldbereich an Volllastkennlinie
% Filtern der nicht anfahrbaren Punkte
for i = 1:length(betr.mot.omega_k_mesh(1,:))
    betr.mot.M_max_mesh(betr.mot.M_max_mesh(:,i)>betr.mot.M_max_vec(i),i) = NaN;
end

% Berechnung der Stroeme (Stator) fuer Punkte im Drehzahl-/Drehmomentgitter
% Minimierung der Stroeme unter Nebenbedingungen (Betriebsbereich Motor)
[ctrl.mot.i_1d_mesh, ctrl.mot.i_1q_mesh] = Optimierung_i(rated, emag, opt, betr.mot.omega_k_mesh, betr.mot.M_max_mesh, 0);

% Berechnung der Stroeme (Rotor) fuer Punkte im Drehzahl-/Drehmomentgitter
ctrl.mot.i_2d_mesh = zeros(size(betr.mot.M_max_mesh));
ctrl.mot.i_2q_mesh = -(emag.L_12./emag.L_22) .* ctrl.mot.i_1q_mesh;

% Berechnen der Spannungen (Stator) fuer Punkte im Drehzahl-/Drehmomentgitter
ctrl.mot.u_1d_mesh = emag.R_1.*ctrl.mot.i_1d_mesh - betr.mot.omega_k_mesh.*emag.L_11.*emag.sigma.*ctrl.mot.i_1q_mesh;
ctrl.mot.u_1q_mesh = emag.R_1.*ctrl.mot.i_1q_mesh + betr.mot.omega_k_mesh.*emag.L_11.*ctrl.mot.i_1d_mesh;

% Berechnen der el. Leistung
ctrl.mot.P_el_mesh = 1.5.*(ctrl.mot.i_1d_mesh.*ctrl.mot.u_1d_mesh + ctrl.mot.i_1q_mesh.*ctrl.mot.u_1q_mesh);

% Berechnung Drehzahlvektoren
ctrl.mot.i_1d_mesh(end,1) = 0;
ctrl.mot.i_1q_mesh(end,1) = 0;
var1 = ctrl.mot.i_1d_mesh(end,:);
var2 = ctrl.mot.i_1q_mesh(end,:);
if(all(ctrl.mot.i_1d_mesh(end,:)<1e-10) && all(ctrl.mot.i_1q_mesh(end,:)<1e-10)) % fuer Berechnung der Winkelgeschwindigkeit bei M = 0 -> i_1d = 0
    ctrl.mot.i_1d_mesh(end,:) = 1e-10;
    ctrl.mot.i_1q_mesh(end,:) = 1e-10;
end
betr.mot.omega_el_mesh = betr.mot.omega_k_mesh - ((emag.R_2./emag.L_22).*(ctrl.mot.i_1q_mesh./ctrl.mot.i_1d_mesh)); % elektrische Winkelgeschwindigkeit
betr.mot.omega_el_mesh(end,1) = 0;
ctrl.mot.i_1d_mesh(end,:) = var1;
ctrl.mot.i_1q_mesh(end,:) = var2;
betr.mot.f_el_mesh = betr.mot.omega_el_mesh ./ (2.*pi);  % elektrische Frequenz
betr.mot.omega_m_mesh = (betr.mot.omega_el_mesh ./ rated.p); % mechanische Winkelgeschwindigkeit
betr.mot.n_m_mesh = (betr.mot.omega_m_mesh ./ (2.*pi)) .* 60; % mechanische Drehzahl
for i = 1:length(betr.mot.n_m_mesh(1,:))
    var = betr.mot.n_m_mesh(~isnan(betr.mot.n_m_mesh(:,i)),i);
    betr.mot.n_m_vec(i) = var(1);
end
clear var1 var2 var

%% B.2) Generatormodell
% #########################################################################
% #   B.2) GENERATORMODELL                                                #
% #########################################################################

if(opt.Generator)
    % Berechnung der Volllastkennlinie (maximales Moment fuer Drehzahlvektor berechnen)
    % Maximierung des Moments unter Nebenbedingungen (Betriebsbereich Generator)
    [ctrl.gen.i_1, ctrl.gen.beta] = Optimierung_M(rated, emag, opt, betr.omega_k_vec, 1);
    betr.gen.M_max_vec = 1.5.*rated.p.*((emag.L_12.^2./emag.L_22).*ctrl.gen.i_1.^2.*0.5.*sin(2.*ctrl.gen.beta));
    betr.gen.M_max_vec = betr.gen.M_max_vec';
    
    % Toleranz für max Moment (Grund: numerische Ungenauigkeit)
    tol = 1e-3;
    
    % Erstellen Momentenvektor
    betr.gen.M_vec = linspace(0, min(betr.gen.M_max_vec)+tol, opt.M_tics);
    clear tol

    % Erstellen Drehzahl-/Drehmomentgitter
    [betr.gen.omega_k_mesh, betr.gen.M_max_mesh] = meshgrid(betr.omega_k_vec, betr.gen.M_vec);
    
    % Anpassung Kennfeldbereich an Volllastkennlinie
    % Filtern der nicht anfahrbaren Punkte
    for i = 1:length(betr.gen.omega_k_mesh(1,:))
        betr.gen.M_max_mesh(betr.gen.M_max_mesh(:,i)<betr.gen.M_max_vec(i),i) = NaN;
    end
    
    % Berechnung der Stroeme (Stator) fuer Punkte im Drehzahl-/Drehmomentgitter
    % Minimierung der Stroeme unter Nebenbedingungen (Betriebsbereich Generator)
    [ctrl.gen.i_1d_mesh, ctrl.gen.i_1q_mesh] = Optimierung_i(rated, emag, opt, betr.gen.omega_k_mesh, betr.gen.M_max_mesh, 1);
    
    % Berechnung der Stroeme (Rotor) fuer Punkte im Drehzahl-/Drehmomentgitter
    ctrl.gen.i_2d_mesh = zeros(size(betr.gen.M_max_mesh));
    ctrl.gen.i_2q_mesh = -(emag.L_12./emag.L_22) .* ctrl.gen.i_1q_mesh;

    % Berechnen der Spannungen (Stator) fuer Punkte im Drehzahl-/Drehmomentgitter
    ctrl.gen.u_1d_mesh = emag.R_1.*ctrl.gen.i_1d_mesh - betr.gen.omega_k_mesh.*emag.L_11.*emag.sigma.*ctrl.gen.i_1q_mesh;
    ctrl.gen.u_1q_mesh = emag.R_1.*ctrl.gen.i_1q_mesh + betr.gen.omega_k_mesh.*emag.L_11.*ctrl.gen.i_1d_mesh;

    % Berechnen der el. Leistung
    ctrl.gen.P_el_mesh = 1.5.*(ctrl.gen.i_1d_mesh.*ctrl.gen.u_1d_mesh + ctrl.gen.i_1q_mesh.*ctrl.gen.u_1q_mesh);

    % Berechnung Drehzahlvektoren
    ctrl.gen.i_1d_mesh(end,1) = 0;
    ctrl.gen.i_1q_mesh(end,1) = 0;
    var1 = ctrl.gen.i_1d_mesh(1,:);
    var2 = ctrl.gen.i_1q_mesh(1,:);
    if(all(ctrl.gen.i_1d_mesh(1,:)<1e-10) && all(ctrl.gen.i_1q_mesh(1,:)<1e-10)) % fuer Berechnung der Winkelgeschwindigkeit bei M = 0 -> i_1d = 0
        ctrl.gen.i_1d_mesh(1,:) = 1e-10;
        ctrl.gen.i_1q_mesh(1,:) = 1e-10;
    end
    betr.gen.omega_el_mesh = betr.gen.omega_k_mesh - ((emag.R_2./emag.L_22).*(ctrl.gen.i_1q_mesh./ctrl.gen.i_1d_mesh)); % elektrische Winkelgeschwindigkeit
    betr.gen.omega_el_mesh(1,1) = 0;
    ctrl.gen.i_1d_mesh(1,:) = var1;
    ctrl.gen.i_1q_mesh(1,:) = var2;
    betr.gen.f_el_mesh = betr.gen.omega_el_mesh ./ (2.*pi);  % elektrische Frequenz
    betr.gen.omega_m_mesh = (betr.gen.omega_el_mesh ./ rated.p); % mechanische Winkelgeschwindigkeit
    betr.gen.n_m_mesh = (betr.gen.omega_m_mesh ./ (2.*pi)) .* 60; % mechanische Drehzahl
    for i = 1:length(betr.gen.n_m_mesh(1,:))
        var = betr.gen.n_m_mesh(~isnan(betr.gen.n_m_mesh(:,i)),i);
        betr.gen.n_m_vec(i) = var(end);
    end
    clear var1 var2 var
    
    % Zusammenfuegen Motor und Generator
    % ctrl.i_1 = [ctrl.mot.i_1; ctrl.gen.i_1(2:end,:)];
    % ctrl.beta = [ctrl.mot.beta; ctrl.gen.beta(2:end,:)];
    % betr.M_max_vec = [betr.mot.M_max_vec betr.gen.M_max_vec(:,2:end)];
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
    % Zusammenfuegen Motor und Generator
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

%% C) Verlustberechnung
% #########################################################################
% #   C) VERLUSTBERECHNUNG                                                #
% #########################################################################

% Berechnung der Wicklungsverluste P_vw [W] (Kupferverluste/Stromwaermeverluste)
% Quelle: [Mueller08, S.438 - Formel 6.3.18]
% die 1.5 kommen aus der Clarke-Transformation, der Strom wird so energierichtig
% (i_d_mesh.^2 + i_q_mesh.^2) sind peak Werte
if(opt.P_vw==1)
    loss.vw.P_1vw_mesh = 1.5 * emag.R_1 * (ctrl.i_1d_mesh.^2 + ctrl.i_1q_mesh.^2);
    loss.vw.P_2vw_mesh = 1.5 * emag.R_2 * (ctrl.i_2d_mesh.^2 + ctrl.i_2q_mesh.^2);
    loss.P_vw_mesh = loss.vw.P_1vw_mesh + loss.vw.P_2vw_mesh;
else
    loss.P_vw_mesh = zeros(size(ctrl.i_1d_mesh));
end

% Berechnung der Ummagnetisierungsverluste P_vu [W]
if(opt.P_vu==1)
    if(strcmp(opt.P_vu_Modell,'Modellansatz Jordan'))
        % Variante 1: Berechnung der Eisenverluste ueber Induktionen und
        % Elektroblechdatenblaetter
        
        % Blechparameter fuer B=1T und f=50Hz berechnen (Stator)
        % Quelle: [Neuschl, S.32-33]
        loss.vu.w_1Fe = opt_Entwurf.Stator_Eisenmaterial.p_vFe(opt_Entwurf.Stator_Eisenmaterial.p_vFe_B_vec==1,:) ./ opt_Entwurf.Stator_Eisenmaterial.p_vFe_f_vec;
        [loss.vu.fitresult_1, ~] = fit(opt_Entwurf.Stator_Eisenmaterial.p_vFe_f_vec', loss.vu.w_1Fe', 'poly1');
        loss.vu.fitresult_data_1 = coeffvalues(loss.vu.fitresult_1);
        loss.vu.sigma_1h_s = loss.vu.fitresult_data_1(2);
        loss.vu.sigma_1w_s = loss.vu.fitresult_data_1(1);
        loss.vu.sigma_1h = loss.vu.sigma_1h_s * 50;
        loss.vu.sigma_1w = loss.vu.sigma_1w_s * 50^2;
        
        % Blechparameter fuer B=1T und f=50Hz berechnen (Rotor)
        % Quelle: [Neuschl, S.32-33]
        loss.vu.w_2Fe = opt_Entwurf.Stator_Eisenmaterial.p_vFe(opt_Entwurf.Rotor_Eisenmaterial.p_vFe_B_vec==1,:) ./ opt_Entwurf.Stator_Eisenmaterial.p_vFe_f_vec;
        [loss.vu.fitresult_2, ~] = fit(opt_Entwurf.Stator_Eisenmaterial.p_vFe_f_vec', loss.vu.w_2Fe', 'poly1');
        loss.vu.fitresult_data_2 = coeffvalues(loss.vu.fitresult_2);
        loss.vu.sigma_2h_s = loss.vu.fitresult_data_2(2);
        loss.vu.sigma_2w_s = loss.vu.fitresult_data_2(1);
        loss.vu.sigma_2h = loss.vu.sigma_2h_s * 50;
        loss.vu.sigma_2w = loss.vu.sigma_2w_s * 50^2;
        
        % Zuschlagsfaktoren
        % Quelle: [Neuschl, S.29 - Formel 4.10]
        loss.vu.c_h = 1;
        loss.vu.c_w = 1;
        
        % Spezifische Verluste in den unterschiedlichen Abschnitten berechnen
        % Quelle: [Neuschl, S.28 - Formel 4.9]
        B_1r = emag.B_1r * ones(size(betr.f_el_mesh));
        B_1z_m = emag.B_1z_m * ones(size(betr.f_el_mesh));
        loss.vu.p_1r = (loss.vu.sigma_1h.*loss.vu.c_h.*(betr.f_el_mesh./50) + loss.vu.sigma_1w.*loss.vu.c_w.*(betr.f_el_mesh./50).^2) .* (B_1r./1).^2;
        loss.vu.p_1z = (loss.vu.sigma_1h.*loss.vu.c_h.*(betr.f_el_mesh./50) + loss.vu.sigma_1w.*loss.vu.c_w.*(betr.f_el_mesh./50).^2) .* (B_1z_m./1).^2;
        
        % Zuschlagsfaktoren
        % Quelle: [Mueller08, S.452 - Tabelle 6.4.2]
        loss.vu.k_u_r = 1.7; % zwischen 1.5 und 1.8
        loss.vu.k_u_z = 1.8; % zwischen 1.7 und 2.5

        % Verluste aus den spezifischen Verlusten berechnen
        loss.vu.P_1r = loss.vu.k_u_r .* loss.vu.p_1r .* opt_Entwurf.Stator_Eisenmaterial.rho_Fe .* geo.Vo_1r;
        loss.vu.P_1z = loss.vu.k_u_z .* loss.vu.p_1z .* opt_Entwurf.Stator_Eisenmaterial.rho_Fe .* geo.Vo_1z;
        
        % Gesamte Eisenverluste
        loss.P_vu_mesh = loss.vu.P_1r + loss.vu.P_1z;
    else
        error('Ungueltige Eingabe bei Variable "opt.P_vu_Modell"');
    end
else
    loss.P_vu_mesh = zeros(size(ctrl.i_1d_mesh));
end

% Berechnung der mechanischen Verluste P_vme [W]
% Quelle: [Mueller08, S.433 - Formel 6.2.1]
if(opt.P_vme==1)
    % Faktoren der Gas- und Lagerreibung k_rb [Ws^2/m^4]
    loss.vme.k_rb = 10; % sollte zwischen 5 und 15 liegen

    % Umfangsgeschwindigkeit Rotor [m/s]
    loss.vme.v_2 = (geo.D_2a ./ 2) .* betr.omega_m_mesh;

    loss.P_vme_mesh = loss.vme.k_rb .* geo.D_2a .* (geo.l_i + 0.8^3 .* 0.6 .* geo.tau_1p) .* loss.vme.v_2.^2;
else
	loss.P_vme_mesh = zeros(size(ctrl.i_1d_mesh));
end

% Berechnung der Zusatzverluste P_zus [W]
% Quelle: [DIN EN 60034-2-1, S.xx]
if(opt.P_vzus==1)
    loss.P_vzus_mesh = abs(ctrl.P_el_mesh) .* (0.025 - 0.005*log(rated.P_N*1e-3));
else
    loss.P_vzus_mesh = zeros(size(ctrl.i_1d_mesh));
end

% Berechnung der Gesamtverluste P_ges [W]
loss.P_vges_mesh = loss.P_vw_mesh + loss.P_vu_mesh + loss.P_vme_mesh + loss.P_vzus_mesh;

%% D) Kennfelderstellung
% #########################################################################
% #   D) KENNFELDERSTELLUNG                                               #
% #########################################################################

% Berechnen der mech. Leistung
betr.mot.P_mech_mesh = betr.mot.M_max_mesh .* betr.mot.omega_m_mesh;

% Berechnung Wirkungsgradkennfelder (Betriebsbereich Motor)
[var1, ~] = size(betr.mot.P_mech_mesh);
eta.mot.eta_vw_mesh_alt = betr.mot.P_mech_mesh ./ ctrl.mot.P_el_mesh;
eta.mot.eta_vw_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vw_mesh(1:var1,1:end));
eta.mot.eta_fe_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vu_mesh(1:var1,1:end));
eta.mot.eta_vme_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vme_mesh(1:var1,1:end));
eta.mot.eta_zus_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vzus_mesh(1:var1,1:end));
eta.mot.eta_ges_mesh = betr.mot.P_mech_mesh ./ (betr.mot.P_mech_mesh + loss.P_vges_mesh(1:var1,1:end));

if(opt.Generator)
    % Berechnen der mech. Leistung
    betr.gen.P_mech_mesh = betr.gen.M_max_mesh .* betr.gen.omega_m_mesh;
    
    % Berechnung Wirkungsgradkennfelder (Betriebsbereich Generator)
    eta.gen.eta_vw_mesh_alt = ctrl.gen.P_el_mesh ./ betr.gen.P_mech_mesh;
    eta.gen.eta_vw_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vw_mesh(var1:end,1:end));
    eta.gen.eta_fe_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vu_mesh(var1:end,1:end));
    eta.gen.eta_vme_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vme_mesh(var1:end,1:end));
    eta.gen.eta_zus_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vzus_mesh(var1:end,1:end));
    eta.gen.eta_ges_mesh = ctrl.gen.P_el_mesh ./ (ctrl.gen.P_el_mesh - loss.P_vges_mesh(var1:end,1:end));
    
    % Zusammenfuegen Motor und Generator
    eta.eta_vw_mesh = [eta.mot.eta_vw_mesh; eta.gen.eta_vw_mesh(2:end,:)];
    eta.eta_fe_mesh = [eta.mot.eta_fe_mesh; eta.gen.eta_fe_mesh(2:end,:)];
    eta.eta_vme_mesh = [eta.mot.eta_vme_mesh; eta.gen.eta_vme_mesh(2:end,:)];
    eta.eta_zus_mesh = [eta.mot.eta_zus_mesh; eta.gen.eta_zus_mesh(2:end,:)];
    eta.eta_ges_mesh = [eta.mot.eta_ges_mesh; eta.gen.eta_ges_mesh(2:end,:)];
else
    % Zusammenfuegen Motor und Generator
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

% Mechanische Drehzahl aequidistant verteilen (fuer Export)
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

% Berechnung erfolgreich
opt.Locked = 1;

% Sortieren der structs
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

% Speichern in Analyse
Analyse.Momentensteuerung = ctrl;
Analyse.Verluste = loss;
Analyse.Betriebsdaten = betr;
Analyse.Wirkungsgrad = eta;

Analyse.Optionen = opt;

end

%% F) Hilfsfunktionen
% #########################################################################
% #   F) HILFSFUNKTIONEN                                                  #
% #########################################################################

%% F.1) Optimierung_M
% #########################################################################
% #   F.1) OPTIMIERUNG_M                                                  #
% #########################################################################

function [i_1, beta] = Optimierung_M(rated, emag, opt, omega_vec, Generator)
%Optimierung Maximiert das Moment unter Nebenbedingungen und Grenzen
%   Notation: i_1 = x(1), beta = x(2)
%   Bound (Motor): 0<=i_1<=i_max,  0<=beta<=pi
%   Linear Inequality Constraint: keine
%   Linear Equality Constraint: keine
%   Nonlinear Constraints: siehe Funktion nonlcon_Optimierung_M

if(Generator)
    % Minimierungsfunktion 
    fun = @(x) -1.5.*rated.p.*((emag.L_12.^2./emag.L_22).*-x(1).^2.*0.5.*sin(2.*x(2)));
    
    % Startwert
    x0 = [0 0];
    % Bound Constraints
    lb = [0 -pi];
    ub = [opt.i_1max 0];
else
    % Minimierungsfunktion 
    fun = @(x) -1.5.*rated.p.*((emag.L_12.^2./emag.L_22).*x(1).^2.*0.5.*sin(2.*x(2)));
    
    % Startwert
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
% Optionen fuer Optimierung
options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point','MaxFunctionEvaluations',3e4,'MaxIterations',2000,...
                            'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-6);

% Speicher allokieren
i_1 = zeros(length(omega_vec(1,:)),1);
beta = zeros(length(omega_vec(1,:)),1);

tic
for j = 1:length(omega_vec(1,:))
    omega_k = omega_vec(1,j);
    % Minimierung mit Startwert [0,0]
    [x_sol,~,exitflag,~] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, ...
    @(x) nonlcon(x, opt.u_1max, emag.L_11, emag.sigma, emag.R_1, omega_k), options);
    if(exitflag~=0 && exitflag~=1 && exitflag~=2)
        % Strom wird inf gesetzt wenn kein Ergebnis gefunden wird
        i_1(j) = inf;
        beta(j) = inf;
    else
        % Zuweisen der Werte wenn Ergebnis gefunden wird
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

%% F.2) Optimierung_i
% #########################################################################
% #   F.2) OPTIMIERUNG_I                                                  #
% #########################################################################

function [i_1dsol, i_1qsol] = Optimierung_i(rated, emag, opt, omega_mesh, M_mesh, Generator)
%Optimierung Minimiert die d- und q-Komponente des Stroms unter
%Nebenbedingungen und Grenzen
%   Notation: i_1d = x(1), i_1q = x(2)
%   Bound (Motor): 0<=i_1d<=i_1max, 0<=i_1q<=i_1max
%   Linear Inequality Constraint: keine
%   Linear Equality Constraint: keine
%   Nonlinear Constraints: siehe Funktion nonlcon_Optimierung_i
%   Aus Gruenden der Rechenzeit werden ab Spalte 2 bzw. Zeile 2 Werte
%   vorheriger Stroeme als Startwerte herangezogen.

% Minimierungsfunktion 
fun = @(x) x(1)^2 + x(2)^2;

if(Generator)
    % Startwert
    x0 = [opt.i_1max -opt.i_1max];
    % Bound Constraints
    lb = [0 -opt.i_1max];
    ub = [opt.i_1max 0];
else
    % Startwert
    x0 = [opt.i_1max opt.i_1max];
    % Bound Constraints
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
% Optionen fuer Optimierung
options = optimoptions(@fmincon,'Display','off','Algorithm','sqp','MaxFunctionEvaluations',3e4,'MaxIterations',400,...
                            'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-6);

% Speicher allokieren
i_1dsol = zeros(size(M_mesh));
i_1qsol = zeros(size(M_mesh));

tic
for j = 1:length(omega_mesh(1,:))
    omega_k = omega_mesh(1,j);
    for k = 1:length(M_mesh(:,1))
        M = M_mesh(k,j);
        if(isnan(M))
            % Aussortieren der nicht anfahrbaren Punkte (Volllastkennlinie)
            i_1dsol(k,j) = NaN;
            i_1qsol(k,j) = NaN;
        else
            if(j>1 && ~isnan(i_1dsol(k,j-1)) && ~isnan(i_1qsol(k,j-1)))
                if(k>1 && ~isnan(i_1dsol(k-1,j)))
                    % Minimierung mit Startwerten aus vorheriger Berechnung
                    [x_sol,~,exitflag,~] = fmincon(fun, [i_1dsol(k-1,j) i_1qsol(k,j-1)], A, b, Aeq, beq, lb, ub, ...
                    @(x) nonlcon(x, M, opt.u_1max, opt.i_1max, emag.L_11, emag.L_12, emag.L_22, emag.sigma, emag.R_1, omega_k, rated.p), options);
                else
                    % Minimierung mit Startwerten aus vorheriger Berechnung
                    [x_sol,~,exitflag,~] = fmincon(fun, [i_1dsol(k,j-1) i_1qsol(k,j-1)], A, b, Aeq, beq, lb, ub, ...
                        @(x) nonlcon(x, M, opt.u_1max, opt.i_1max, emag.L_11, emag.L_12, emag.L_22, emag.sigma, emag.R_1, omega_k, rated.p), options);
                end
            else
                % Minimierung mit Startwert [0,0]
                [x_sol,~,exitflag,~] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, ...
                @(x) nonlcon(x, M, opt.u_1max, opt.i_1max, emag.L_11, emag.L_12, emag.L_22, emag.sigma, emag.R_1, omega_k, rated.p), options);
            end
            % Werte < tol werden zu 0 gesetzt (numerische Ungenauigkeit)
            if(abs(x_sol(1))<1e-2)
                x_sol(1)=0;
            end
            if(abs(x_sol(2))<1e-2)
                x_sol(2)=0;
            end
            if(exitflag~=0 && exitflag~=1 && exitflag~=2)
                % Minimierung mit Startwert [0,0]
                [x_sol,~,exitflag,~] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, ...
                @(x) nonlcon(x, M, opt.u_1max, opt.i_1max, emag.L_11, emag.L_12, emag.L_22, emag.sigma, emag.R_1, omega_k, rated.p), options);
                if(exitflag~=0 && exitflag~=1 && exitflag~=2)
                    % Strom wird inf gesetzt wenn kein Ergebnis gefunden wird
                    i_1dsol(k,j) = inf;
                    i_1qsol(k,j) = inf;
                    warning(['fmincon exitflag:', num2str(exitflag), ', check function Optimierung_i'])
                else
                    % Zuweisen der Werte wenn Ergebnis gefunden wird
                    i_1dsol(k,j) = x_sol(1);
                    i_1qsol(k,j) = x_sol(2);
                end
            else
                % Zuweisen der Werte wenn Ergebnis gefunden wird
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
