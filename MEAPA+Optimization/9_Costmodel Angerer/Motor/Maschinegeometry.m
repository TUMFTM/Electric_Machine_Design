function [L_M, ...% Length with housing
            D, ...% armature diameter
            MV, ...% active volume
            D_a, ...% outer diameter of stator
            l_fe, ...% laminated core length
            A_wS, ...% winding cross section stator 
            A_wL, ...% winding cross section rotor
            l_WK, ...% length of winding head
            D_rot, ...% Rotor diameter
             D_aussen_EM, ...% Housing outer diameter
            D_i,...% Rotor inner diameter
             P_n, ...& nominal power 
            V_Gehaeuse, ... %Volume of the housing
            A_nuten_Stator, ...% Stator slot Cross-sectional area
            A_nuten_Rotor, ...% Rotor slot Cross-sectional area
            A_Luftspalt_und_Magnete,...% Cross-sectional area of air gap and magnets
            A_Welle,... %Cross-sectional area of rotor shaft
            N_S, ... % Number of stator slots
            N_L]=... % Number of rotor slots
    Maschinegeometry(M_n, n_n, Maschinentyp, U_linetoline, Maschinendaten)

p=Maschinendaten.Entwurf.Bemessungswerte.p;

[mu_0, p, U_n, m, P_n, Material_Stator, Material_Laeufer, Schaltung, Bauweise, Bauweise_Rotor, Kuehlung, f_n, ...
            w, U_s, E_h, B_m_max, S_S_max, S_L_max, phi_n, phi_fe, B_zS_max, alpha_p_ang, ...
            B_zL_max, B_rS_max, B_rL_max, plothandle, ...
            S_R, lambda, lambdaschalter, eta, etaschalter, ...
            cphi, cphischalter, C_m, cmschalter, ...
            n_s]=...
            Initialisation_EM(M_n, n_n, p, Maschinentyp, U_linetoline, Maschinendaten);

ws = warning('off','all');

%% Rough design with machine equation

% Required active electrical power P_w1 in kW
P_el=Maschinendaten.Entwurf.EMAG.P_el_N;
% Ermittlung des Nennstroms I_n in A
I_n=Maschinendaten.Entwurf.EMAG.I_N;
% Determination of the rated current I_n1 per string (after the circuit) in A
I_n1=Maschinendaten.Entwurf.EMAG.P_s*1000/(U_s*m);
% String current I_s (in one phase of the network) in A
I_s=Maschinendaten.Entwurf.EMAG.P_s*1000/(m*U_n);

    P_i = Maschinendaten.Entwurf.Bemessungswerte.P_N*E_h/(U_s*cphi*Maschinendaten.Entwurf.EMAG.eta_N/100);
    P_n = Maschinendaten.Entwurf.Bemessungswerte.P_N;
     
%    [D,delta,T_p,l_i,l_fe,l ] =  hauptabmessungen(P_n,p,C_m,n_s,Kuehlung,Maschinendaten );
D = Maschinendaten.Entwurf.Geometrie.D_2a;
delta = Maschinendaten.Entwurf.Geometrie.delta;
T_p = Maschinendaten.Entwurf.Geometrie.tau_1p  ;
l_i = Maschinendaten.Entwurf.Geometrie.l_i;
l_fe = Maschinendaten.Entwurf.Geometrie.l_Fe;
l = Maschinendaten.Entwurf.Geometrie.l;

%% winding design
E_h_r = E_h;
a = Maschinendaten.Entwurf.Wicklung.a_1;

w_s_opt = Maschinendaten.Entwurf.EMAG.w_1Str_opt;
q_S = Maschinendaten.Entwurf.Wicklung.q_1;
z_nS = Maschinendaten.Entwurf.Wicklung.z_1n;
w_s = Maschinendaten.Entwurf.EMAG.w_1Str;
N_S = Maschinendaten.Entwurf.Wicklung.N_1;
T_nS = Maschinendaten.Entwurf.Wicklung.tau_1n;
w_ges = w_s*m*a;
z_S = w_ges*2;
q_S_alt =q_S_schleife;
w_s_fehler = Maschinendaten.Entwurf.Wicklung.err_1;

Z = q_S;
xi_1_tat = Maschinendaten.Entwurf.Wicklung.xi_1p;
B_max = B_m_max;
B_m_tat = B_max;
B_m_ang = alpha_p_ang * B_m_tat;
Phi_h = Maschinendaten.Entwurf.EMAG.Phi_h;
Phi_delta = Maschinendaten.Entwurf.EMAG.Phi_delta;
Phi_delta_tat = E_h*sqrt(2)/(w_s*w*xi_1_tat);
if Maschinendaten.Entwurf.Wicklung.Wicklungstyp_1  == '1SGL'
    Wicklungsart ='Einschichtwicklung';
else
   Wicklungsart ='gesehnte Zweischichtwicklung';
end

I_s_zw = Maschinendaten.Entwurf.EMAG.I_1zw;
A_S = Maschinendaten.Entwurf.EMAG.A;
A_L = A_S*cphi;
A_1S = Maschinendaten.Entwurf.Wicklung.A_1L;
A_nS = Maschinendaten.Entwurf.Geometrie.Nut_1.A_1n;
A_wS = A_1S*N_S*z_nS;
A_wL = A_wS*(S_S_max/S_L_max)*cphi;

b_zS = Maschinendaten.Entwurf.Geometrie.Nut_1.b_1z_m;
h_nS = Maschinendaten.Entwurf.Geometrie.Nut_1.h_1n_max;
b_nS = Maschinendaten.Entwurf.Geometrie.Nut_1.b_1n_o;
b_nS_m = Maschinendaten.Entwurf.Geometrie.Nut_1.b_1n_m;

if LDS_values.type == 'ASM'
N_L = Maschinendaten.Entwurf.Wicklung.N_2;
T_nL = Maschinendaten.Entwurf.Wicklung.tau_2n;
A_nL = Maschinendaten.Entwurf.Geometrie.Nut_2.A_2n;
b_zL = Maschinendaten.Entwurf.Geometrie.Nut_2.b_2z_m;
h_nL = Maschinendaten.Entwurf.Geometrie.Nut_2.h_2n_max;
b_nL = Maschinendaten.Entwurf.Geometrie.Nut_2.b_2n_o;
q_L = Maschinendaten.Entwurf.Wicklung.q_2;
h_nL_max = h_nL;
else
end

Phi_rS_max = Maschinendaten.Entwurf.EMAG.Phi_1r_max;
h_rS = Maschinendaten.Entwurf.Geometrie.Nut_1.h_1r_max;
Phi_rL_max = Maschinendaten.Entwurf.EMAG.Phi_2r_max;
h_rL = Maschinendaten.Entwurf.Geometrie.Nut_2.h_2r_max;
D_a = Maschinendaten.Entwurf.Geometrie.D_1a_max;
A_rS = Maschinendaten.Entwurf.Geometrie.A_1r* 1000000; 
A_rL = Maschinendaten.Entwurf.Geometrie.A_2r* 1000000;
l_i_tat = l_i;
T_nL_tat = Maschinendaten.Entwurf.Wicklung.tau_2n;
b_sS = Maschinendaten.Entwurf.Geometrie.Nut_1.b_1ns;
b_sL = Maschinendaten.Entwurf.Geometrie.Nut_2.b_2ns;
gamma_S = 1/(1+5*delta/(b_sS*1000));                          
gamma_L = 1/(1+5*delta/(b_sL*1000));                          
k_cS = T_nS/(T_nS-gamma_S*b_sS);
k_cL = T_nL_tat/(T_nL_tat-gamma_L*b_sL);
k_c = k_cS*k_cL;
V_delta_ang = Maschinendaten.Entwurf.EMAG.V_delta;
V_delta_tat = V_delta_ang;
B_zS_abg = Maschinendaten.Entwurf.EMAG.B_1z_m;
B_zL_abg = Maschinendaten.Entwurf.EMAG.B_2z_m;
D_i = Maschinendaten.Entwurf.Geometrie.D_1i;


t_gehaeuse=0.01; %in m
L_gehaeuseueberhang=0.0459; %in m
l_WK=0.03; %in m
L_M=l_fe+2*l_WK+L_gehaeuseueberhang;
MV =(D_a^2*pi/4)*L_M;   
D_rot=D;
t_luft=0.005;
D_aussen_EM=D_a+2*(t_gehaeuse+t_luft);
V_Gehaeuse=((D_aussen_EM^2-(D_aussen_EM-2*t_gehaeuse)^2)*pi/4)*(l_fe+2*l_WK)+D_aussen_EM^2*pi/4*2*t_gehaeuse;
A_nuten_Stator=A_nS*N_S/1000000;
A_nuten_Rotor=A_nL*N_L/1000000;
A_Luftspalt_und_Magnete=((D+0.0085)^2*pi/4)-(D^2*pi/4);

alpha_p_tat = alpha_p_ang;

% Estimated voltage drop in the teeth in A
V_zS_abg=Maschinendaten.Entwurf.EMAG.V_1z;
V_zL_abg=Maschinendaten.Entwurf.EMAG.V_2z;
% Saturation factor k in 1
k=(V_zS_abg+V_zL_abg)/V_delta_ang;

B_zS_tat=1.03*B_max_tat*T_nS*l_i_tat/(l_fe*phi_fe*b_zS);
B_zL_tat=B_max_tat*T_nL_tat*l_i_tat/(l_fe*phi_fe*b_zL);

[d_welle_a, d_welle_i]=Rotor_Design(M_n, D_i);

A_Welle=(d_welle_a^2-d_welle_i^2)*pi/4;


P_n_VA=P_n;

