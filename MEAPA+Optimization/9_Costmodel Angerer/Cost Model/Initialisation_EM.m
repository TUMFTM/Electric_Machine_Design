function [mu_0, p, U_n, m, P_n, Material_Stator, Material_Laeufer, Schaltung, Bauweise, Bauweise_Rotor, Kuehlung, f_n, ...
            w, U_s, E_h, B_m_max, S_S_max, S_L_max, phi_n, phi_fe, B_zS_max, alpha_p, ...
            B_zL_max, B_rS_max, B_rL_max, plothandle, ...
            S_R, lambda, lambdaschalter, eta, etaschalter, ...
            cphi, cphischalter, C_m, cmschalter, ...
            n_s]=...
            Initialisation_EM(M_n, n_n, p, Maschinentyp, U_linetoline, Maschinendaten) 
        
        

mu_0=Maschinendaten.Entwurf.Konstanten.mu_0;
%% Reading in the motor data
U_n=U_linetoline;         % in V         
m=Maschinendaten.Entwurf.Bemessungswerte.m;                             
P_n=Maschinendaten.Entwurf.Bemessungswerte.P_N;    
%% Material data
if Maschinendaten.Entwurf.Optionen.Stator_Leitermaterial.Bezeichnung == 'Kupfer'
    Material_Stator='Kupfer';
    Material_Laeufer='Kupfer' ;
else
    disp('wrong');
end
%% circuit data
if Maschinendaten.Entwurf.Optionen.Schaltung == 'Stern'
    Schaltung='Stern';
else
    disp('wrong');
end
%% Topology
Bauweise='geschlossen';
Bauweise_Rotor='Schenkelpol';
%% Cooling type
if strcmp(Maschinendaten.Entwurf.Optionen.Kuehlungsart,'Wasser (direkt)') || strcmp(Maschinendaten.Entwurf.Optionen.Kuehlungsart,'Innen- oder Kreislaufkuehlung') 
Kuehlung='Fluessig';
else
    Kuehlung='Luft';
    disp('Falsche Kühlungsart');
end
%% Angular velocity of the supply in rad/s
f_n= Maschinendaten.Entwurf.Bemessungswerte.f_N;
w=2*pi*f_n;
n_s=f_n/p*60;                    % [1/min]


U_s = Maschinendaten.Entwurf.EMAG.U_1Str;

%% Determination of the induced voltage E_h in V
E_h = Maschinendaten.Entwurf.EMAG.E_h;

%% Read variables for winding design, slot shapes and magnetic circuit

if Maschinentyp == 'ASM'
    B_m_max = Maschinendaten.Entwurf.EMAG.B_m;
else
    B_m_max = Maschinendaten.Entwurf.Richtwerte.B_m;
end

% Ring current density
S_R = Maschinendaten.Entwurf.Wicklung.S_2r; 

% maximum current density in the stator in A/mm^2
S_S_max=Maschinendaten.Entwurf.Wicklung.S_1;

% maximum current density S_l in the rotor in A/mm^2
S_L_max = Maschinendaten.Entwurf.Wicklung.S_2;

% Groove filling factor
phi_n=Maschinendaten.Entwurf.Richtwerte.phi_1n;
phi_fe=Maschinendaten.Entwurf.Richtwerte.phi_1Fe;

% if strcmp(Maschinentyp,'PMSM')==1
B_zS_max=Maschinendaten.Entwurf.Richtwerte.B_1z_max;                                  
B_zL_max=Maschinendaten.Entwurf.Richtwerte.B_2z_max;                                 
B_rS_max=Maschinendaten.Entwurf.Richtwerte.B_1r_max;                                       
B_rL_max=Maschinendaten.Entwurf.Richtwerte.B_2r_max;                                       

alpha_p=Maschinendaten.Entwurf.EMAG.alpha_p;

plothandle = 'NEIN';

%% Switch for presetting certain values
lambda = Maschinendaten.Entwurf.Richtwerte.lambda; 
lambdaschalter = 1;
eta=Maschinendaten.Entwurf.EMAG.eta_N;          
etaschalter = 1; 
cphi=Maschinendaten.Entwurf.EMAG.cos_phi_N;            
cphischalter = 1; 
if vehicle.Maschinendaten.type == 'PMSM'
    C_m = Maschinendaten.Entwurf.Geometrie.misc.C_s; 
else
    C_m = Maschinendaten.Entwurf.Geometrie.misc.C_mech;
end
cmschalter = 1;





