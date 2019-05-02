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

function [Verluste] = Verluste_PMSM(Maschinendaten, i_d_mesh, i_q_mesh, omega_k_mesh)
% This function calculates the losses of a PMSM, based on the selected losses.
% The iron losses are based on the losses the sheet metal manufacturer specify 
% on their data sheets.

%% Parameter re-storage for easier use
% Rated sizes
prim = Maschinendaten.Bemessungsgroessen.Primaerparameter;
sek = Maschinendaten.Bemessungsgroessen.Sekundaerparameter;
% Approx. values
richt = Maschinendaten.Richtwerte;
% Options
opt = Maschinendaten.Optionen;
% Design
ent = Maschinendaten.Entwurf;
% Losses
sel = Maschinendaten.Verluste.Auswahl;

%% Calculation of winding losses P_vw [W] (copper losses/current heat losses)
% Literature: [Mueller08, p.438 - Formula 6.3.18]
% The 1.5 comes from the Clarke-transformation.
% (i_d_mesh.^2 + i_q_mesh.^2) represent peak values
if(sel.Stromwaermeverluste==1)
    P_vw_mesh = 1.5 * ent.R_s * (i_d_mesh.^2 + i_q_mesh.^2);
else
    P_vw_mesh = zeros(size(i_d_mesh));
end

%% Calculation of the iron losses P_fe [W]
if(sel.Eisenverluste==1)
   
    % Literature: Lecture notes Drive control for electric vehicles
    P_fe_mesh = 1.5.*omega_k_mesh.^(1.5).*((ent.L_d.*i_d_mesh + ent.psi_PM).^2 + (ent.L_q.*i_q_mesh).^2);
else
    P_fe_mesh = zeros(size(i_d_mesh));
end

%% Calculation of the mechanical losses P_vme [W]
% Literature: [Mueller08, p.433 - Formula 6.2.1]
if(sel.mechanische_Verluste==1)
    % Factors of air and bearing friction k_rb [Ws^2/m^4]
    k_rb = 10; % should be between 15 and 5

    % Circumferential speed Rotor [m/s]
    v_2 = ent.d_a / 2 * (omega_k_mesh ./ prim.p);

    P_vme_mesh = k_rb * ent.d_a * (ent.l_i + 0.8^3 * 0.6 * ent.tau_p) * v_2.^2;
else
	P_vme_mesh = zeros(size(i_d_mesh));
end

%% Calculation of mechanical losses P_vme alternatively [W]
% Literature: [Finken11, p.71 - Formula 5.44]
% Note: equation above  is more accurate
% P_vme_mesh = 0.01 * prim.P_N * 1e3 * ((omega_k_mesh ./ (2 .* pi) ./ prim.p .* 60) ./ prim.n_N);

%% Calculation of the additional losses P_zus [W]
% Literature: [DIN EN 60034-2-1]

if(sel.Zusatzverluste==1)
    P_zus_mesh = zeros(size(i_d_mesh));
else
    P_zus_mesh = zeros(size(i_d_mesh));
end

%% Calculation of the total losses P_ges [W]
P_vges_mesh = P_vw_mesh + P_fe_mesh + P_vme_mesh + P_zus_mesh;

%% Save the calculated parameters to machine data struct
Verluste = struct('P_vw_mesh',P_vw_mesh,'P_fe_mesh',P_fe_mesh, ...
    'P_vme_mesh',P_vme_mesh,'P_zus_mesh',P_zus_mesh,'P_vges_mesh',P_vges_mesh,'Auswahl',sel);

end