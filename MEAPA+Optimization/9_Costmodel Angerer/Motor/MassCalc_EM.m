function [m_gesamt, m_alu, m_elektroblech, m_kupferdraht, m_magnet, m_rotorkaefig, m_stahl]=MassCalc_EM(D, D_i, MV, A_wS, A_wL, l_WK, D_a, l_fe, D_rot, A_nuten_Stator, A_nuten_Rotor, A_Luftspalt, A_Magnete, A_Welle, typ_EM, V_Gehaeuse, Maschinendaten)
%% Material data

Di_EB=7600; %density electrical sheet in kg/m^3
Di_KD=8920; %density copper wire in kg/m^3
Di_AL=2700; %density aluminum in kg/m^3
Di_ST=7700; %density steel alloy in kg/m^3
Di_MG=7400; %density magnetic material in kg/m^3
Stapelfaktor=0.96;  %stacking factor minimum 0.92 according to telephone conversation with Sonja Flügel from the Herzog Chair; according to DIN 46400, the S. according to Tab. l is between 0.95 and 0.97


% Stator and rotor material depending on magnet arrangement
if strcmp(typ_EM,'ASM') && strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'Kaefiglaeufer')
m_elektroblech= (((D_a^2)/4*pi)... % Stator outside diameter used, for total machine in m
                    -A_nuten_Stator * 1e-6... % stator slots subtracted, in mm^2, therefore * 1e-6, now in m^2
                    -A_Luftspalt... % Cross-sectional area of air gap and magnets in m, actually only air gap
                    -((D_i^2)/4*pi)... % Axis is taken out in m^2
                    -(A_wL/1000000))...     % Rotor solid material minus part of the cage rotor by active length, skewing not considered here!
                    *l_fe*Di_EB*Stapelfaktor; %*0.7; in kg
                
elseif strcmp(typ_EM,'ASM') && strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'Schleifringlaeufer')
m_elektroblech= (((D_a^2)/4*pi)... % Stator outside diameter used, for total machine in m
                    -A_nuten_Stator * 1e-6... % stator slots subtracted, in mm^2, therefore * 1e-6, now in m^2
                    -A_nuten_Rotor * 1e-6...  % rotor slots subtracted, in mm^2, therefore * 1e-6, now in m^2
                    -A_Luftspalt... % Cross-sectional area of air gap and magnets in m, actually only air gap
                    -((D_i^2)/4*pi))... % Axis is taken out in m^2
                    *l_fe*Di_EB*Stapelfaktor; %*0.7 in kg
                
elseif strcmp(typ_EM,'PMSM') && strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'SPMSM') %Bei SPMSM sind Magnete in Luftspalt abgebildet
m_elektroblech= (((D_a^2)/4*pi)... % Stator outside diameter used, for total machine in m
                    -A_nuten_Stator * 1e-6... % stator slots subtracted, in mm^2, therefore * 1e-6, now in m^2
                    -A_Luftspalt... % Cross-sectional area of air gap and magnets in m, actually only air gap
                    -((D_i^2)/4*pi))... % Axis is taken out in m^2
                    *l_fe*Di_EB*Stapelfaktor; %*0.7 in kg

elseif strcmp(typ_EM,'PMSM') && ((strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'IPMSM (tangential)')) || (strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'IPMSM (V-Form)')) || (strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'IPMSM (eingelassen)')))
m_elektroblech= (((D_a^2)/4*pi)... % Stator outside diameter used, for total machine in m
                    -A_nuten_Stator * 1e-6... % stator slots subtracted, in mm^2, therefore * 1e-6, now in m^2
                    -A_Luftspalt... % Cross-sectional area of air gap and magnets in m, actually only air gap
                    -((D_i^2)/4*pi)... % Axis is taken out in m^2
                    -A_Magnete)... %in m^2
                    *l_fe*Di_EB*Stapelfaktor; %*0.7 in kg 
else
    printf('Error in mass calculation for typ_EM:', typ_EM, ' and Maschienenausfuehrung:', Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung);
end
    
% Housing   
    m_alu=V_Gehaeuse*Di_AL;   %in kg

% shaft
    m_stahl=A_Welle*2*l_fe*Di_ST;   %in kg
 
% Copper wire,rotor cage and magnets!
if strcmp(typ_EM,'PMSM')
    m_magnet=(A_Magnete)*l_fe*Di_MG; 
    m_kupferdraht=(A_wS/1000000)*(2*l_WK+l_fe)*Di_KD;   %in kg
    m_rotorkaefig=0;
elseif strcmp(typ_EM,'ASM') && strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'Kaefiglaeufer')
    m_rotorkaefig=((A_wL*l_fe/1000000)+(2*D^2*pi/4*0.01))*Di_AL;
    m_kupferdraht=((A_wS/1000000)*(2*l_WK+l_fe))*Di_KD;   %in kg
    m_magnet=0;
elseif strcmp(typ_EM,'ASM') && strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'Schleifringlaeufer')
    m_kupferdraht=((A_wS/1000000)*(2*l_WK+l_fe)+(A_wL/1000000)*(2*l_WK+l_fe))*Di_KD;   %in kg
    m_magnet=0;
end

% Factor 0.9: In reality, the rotor is not fitted over the entire surface 
% because certain set points for robot gripper arms or handling devices are necessary. 
% handling devices are necessary. This factor was used to 
% was approximately taken into account. In logic, the 10 % corresponds to the 
% material loss of the purchased masses, since my assumption was, 
% that material is always lost. 
% But not needed here, because Tool outputs correct value

m_gesamt=(m_alu+m_elektroblech+m_kupferdraht+m_magnet+m_rotorkaefig+m_stahl); % in kg



