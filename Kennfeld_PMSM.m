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

function [Maschinendaten, handles] = Kennfeld_PMSM(Maschinendaten)
% The Map_PMSM function calculates efficiency diagrams for the previously
% calculated PMSM. The currents and voltages
% are calculated using a motor model. Then the individual
% Note: There is no thermal check.
% The generator operation is not implemented.
% V/A: linear & stationary consideration:
% - no saturation
% - no cross coupling of inductances
% - Inductances are not current-dependent
% - no dynamics
% - interlinked flow of PM constant

%% Parameter re-storage for easier use
% Rated values
prim = Maschinendaten.Bemessungsgroessen.Primaerparameter;
sek = Maschinendaten.Bemessungsgroessen.Sekundaerparameter;
% Approx. values
richt = Maschinendaten.Richtwerte;
% Options
opt = Maschinendaten.Optionen;
% Design
ent = Maschinendaten.Entwurf;

%% Start of the GUI
handles = GUI_Kennfeld_PMSM(Maschinendaten);
if(~isstruct(handles))
    error('Process aborted by user!');
else
    % Write struct in machine data
    Maschinendaten.Regelgroessen = handles.Regelgroessen;
    reg = Maschinendaten.Regelgroessen;
    Maschinendaten.Verluste.Auswahl = handles.Verluste;
    Maschinendaten.Optionen.axes_plot = handles.axes_plot;
    names = [fieldnames(Maschinendaten.Optionen); fieldnames(handles.Optionen)];
    Maschinendaten.Optionen = cell2struct([struct2cell(Maschinendaten.Optionen); struct2cell(handles.Optionen)], names, 1);
    opt = Maschinendaten.Optionen;
end

%% Preparation of the map area
% Check max speed
if(opt.n_max_Value==1)
    omega_k_max = (opt.n_max/60)*(2*pi)*prim.p;
else
    omega_k_max = opt.n_max*(prim.n_N/60)*(2*pi)*prim.p;
end

i_f = ent.psi_PM / ent.L_d;
if(i_f>=reg.i_max)
    omega_k_max2 = reg.u_max / (ent.psi_PM - ent.L_d*reg.i_max) - 5;
    if(omega_k_max>omega_k_max2)
        omega_k_max = omega_k_max2;
        warning('max. speed limited by high interlinked flux of PM');
    end
end

% Create speed vector
omega_k_vec = linspace(0, omega_k_max, opt.tics_omega);

%% Calculating the currents for operating points in the speed/torque grid
% This function calculates the full load characteristic and currents of the PMSM.
% More detailed information on the calculations can be found in the function.
[i_d_mesh, i_q_mesh, M_max_vec, M_max_mesh, omega_k_mesh] = Motormodell(prim, ent, reg, omega_k_vec, opt.tics_M);

%% Calculating the voltages for operating points in the speed/torque grid
u_d_mesh = ent.R_s.*i_d_mesh - omega_k_mesh.*ent.L_q.*i_q_mesh;
u_q_mesh = ent.R_s.*i_q_mesh + omega_k_mesh.*(ent.L_d.*i_d_mesh + ent.psi_PM);

%% Calculation of mech. power
P_mech_mesh = M_max_mesh .* (omega_k_mesh/prim.p);

%% Calculation of the el. power
P_el_mesh = 1.5.*(i_d_mesh.*u_d_mesh + i_q_mesh.*u_q_mesh);

%% Loss calculation PMSM
% This function calculates the losses of the PMSM. More information about 
% the calculations can be found in the function.
[Verluste] = Verluste_PMSM(Maschinendaten, i_d_mesh, i_q_mesh, omega_k_mesh);

%% Calculation of efficiency characteristic maps
eta_mesh = P_mech_mesh ./ P_el_mesh;
eta_vw_mesh = P_mech_mesh ./ (P_mech_mesh + Verluste.P_vw_mesh);
eta_fe_mesh = P_mech_mesh ./ (P_mech_mesh + Verluste.P_fe_mesh);
eta_vme_mesh = P_mech_mesh ./ (P_mech_mesh + Verluste.P_vme_mesh);
eta_zus_mesh = P_mech_mesh ./ (P_mech_mesh + Verluste.P_zus_mesh);
eta_ges_mesh = P_mech_mesh ./ (P_mech_mesh + Verluste.P_vges_mesh);

%% Save the calculated parameters to machine data struct
Maschinendaten = rmfield(Maschinendaten,'Verluste');
Maschinendaten.Kennfeld = struct('omega_k_vec',omega_k_vec,'omega_k_mesh',omega_k_mesh,...
    'M_max_vec',M_max_vec,'M_max_mesh',M_max_mesh,'P_mech_mesh',P_mech_mesh,...
    'P_el_mesh',P_el_mesh,'i_d_mesh',i_d_mesh,'i_q_mesh',i_q_mesh,'u_d_mesh',u_d_mesh,...
    'u_q_mesh',u_q_mesh,'Verluste',Verluste,'eta_vw_mesh',eta_vw_mesh,...
    'eta_fe_mesh',eta_fe_mesh,'eta_vme_mesh',eta_vme_mesh,'eta_zus_mesh',eta_zus_mesh,...
    'eta_ges_mesh',eta_ges_mesh);

%% Feedback GUI
handles.pushbutton_start.FontSize = 11;
handles.pushbutton_start.String = '<html><center>Kennfeldberechnung erfolgreich abgeschlossen<br>Ergebnisse in Maschinendaten gespeichert';
handles.popupmenu_Auswahl_Plot.Enable = 'on';

var = {'Please choose ...'};
if(handles.Verluste.Stromwaermeverluste==1 | handles.Verluste.Eisenverluste==1 |...
       handles.Verluste.mechanische_Verluste==1 |  handles.Verluste.Zusatzverluste==1)
   var = [var, 'Wirkungsgrad gesamt', 'Verluste gesamt'];
   if(handles.Verluste.Stromwaermeverluste==1)
        var = [var, 'Wirkungsgrad Stromwaermeverluste','Stromwaermeverluste'];
    end
    if(handles.Verluste.Eisenverluste==1)
        var = [var, 'Wirkungsgrad Eisenverluste', 'Eisenverluste'];
    end
    if(handles.Verluste.mechanische_Verluste==1)
        var = [var, 'Wirkungsgrad mechanische Verluste', 'mechanische Verluste'];
    end
    if(handles.Verluste.Zusatzverluste==1)
        var = [var, 'Wirkungsgrad Zusatzverluste', 'Zusatzverluste'];
    end
    var = [var, 'i_d', 'i_q', 'u_d', 'u_q'];
else
    var = [var,'i_d', 'i_q', 'u_d', 'u_q'];
end
handles.popupmenu_Auswahl_Plot.String = var;

var = msgbox({'Map creation successfully completed'},'Success','help','modal');
set(var,'Position',[1081 673 259 80])
waitfor(var);

%% Save struct Maschinendaten
save(['Ergebnisse/',Maschinendaten.Optionen.folder_id,'/Maschinendaten_',Maschinendaten.Optionen.file_id,'.mat'],'Maschinendaten');

end