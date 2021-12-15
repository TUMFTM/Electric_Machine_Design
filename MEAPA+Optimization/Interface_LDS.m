% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Modell for the design and analysis of PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (svenja.kalt@tum.de), 
%           Jonathan Erhard 
% -------------------------------------------------------------------------

function motor = Interface_LDS(Maschinendaten)
% Die Funktion Interface_LDS konvertiert die Daten aus der Kennfeldermittlung der PMSM in die Form für die LDS Berechnung von Krapf/Koch.
 
%% Convert data to correct form
% Design sizes
motor.type = Maschinendaten.Entwurf.Optionen.Maschinentyp;
motor.power = Maschinendaten.Entwurf.Bemessungswerte.P_N/1000;                                                                % Rated power in kW
motor.n_n = Maschinendaten.Entwurf.Bemessungswerte.n_N;                                                                       % Rated rotational speed in U/min
motor.n_max = Maschinendaten.Analyse.Optionen.n_max;                                                                                                                      %Max. Drehzahl in U/min
motor.U_n = Maschinendaten.Entwurf.Bemessungswerte.U_N;                                                                       % rotational voltage in V

% Speed axis - X-axis map
motor.eff_n_axis = transpose(Maschinendaten.Analyse.Betriebsdaten.mot.n_m_vec);

% Torque - Y-Axis
Maschinendaten.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(Maschinendaten.Analyse.Betriebsdaten.mot.M_vec));
motor.eff_T_axis = [(fliplr(transpose(Maschinendaten.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(Maschinendaten.Analyse.Betriebsdaten.mot.M_vec)];        %Original-Variable
%motor.eff_T_axis = [flipud(transpose(Maschinendaten.Analyse.Betriebsdaten.mot.M_vec));(fliplr(transpose(Maschinendaten.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1)];        %Original-Variable

%Diagram                                                            
    % Manipulate efficiency map
    Maschinendaten.Analyse.Wirkungsgrad.eta_ges_mesh=flipud(Maschinendaten.Analyse.Wirkungsgrad.eta_ges_mesh);
    Maschinendaten.Analyse.Wirkungsgrad.eta_ges_mesh(Maschinendaten.Analyse.Wirkungsgrad.eta_ges_mesh<0.1) = 0.1; %set small efficiency values to 0.1
    Maschinendaten.Analyse.Wirkungsgrad.eta_ges_mesh(1,:) = 0; %first colomn is zero
    Maschinendaten.Analyse.Wirkungsgrad.eta_ges_mesh(:,1) = 0.1; %first row is zero
    % Create one "larger" eff. map that includes acceleration and recuperation
    motor.eff_recu = 1./Maschinendaten.Analyse.Wirkungsgrad.eta_ges_mesh; %recuperation efficiency
    motor.eff_recu(motor.eff_recu==inf)=100;
    motor.eff_recu = motor.eff_recu(2:end,:);
    motor.eff = [flipud(motor.eff_recu); Maschinendaten.Analyse.Wirkungsgrad.eta_ges_mesh];

% full load characteristic
%motor.T_max = transpose (Maschinendaten.Analyse.Betriebsdaten.mot.M_max_vec);                                                            %noch auf 2000 interpolieren?
Maschinendaten.Analyse.Betriebsdaten.mot.M_max_vec=transpose(Maschinendaten.Analyse.Betriebsdaten.mot.M_max_vec);
motor.T_max = transpose(Maschinendaten.Analyse.Betriebsdaten.mot.M_max_vec);  

% Speed axis scaled
motor.T_max_n_axis = transpose(Maschinendaten.Analyse.Betriebsdaten.mot.M_vec);                                                                          %noch auf 2000 interpolieren?

Verlustleistung = Maschinendaten.Analyse.Verluste.P_vges_mesh;
Graphtitle = 'Verlustleistungskennfeld';
folder_id = [Maschinendaten.Entwurf.Optionen.folder_id]; %file_id, '_data'];
file_id = [Maschinendaten.Entwurf.Optionen.file_id];
% save (['3_Ergebnisse/',folder_id,'/','2_Analyse/', motor.type,'_',num2str(motor.power),'_',num2str(motor.n_n),'_',num2str(motor.n_max),'_',num2str(motor.U_n),'.mat'],'motor');
% save (['3_Ergebnisse/',folder_id,'/','2_Analyse/', Graphtitle, '_' ,folder_id,'.mat'],'Verlustleistung');
% save (['LDS/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', motor.type,'_',num2str(motor.power),'_',num2str(motor.n_n),'_',num2str(motor.n_max),'_',num2str(motor.U_n),'.mat'],'motor');

end