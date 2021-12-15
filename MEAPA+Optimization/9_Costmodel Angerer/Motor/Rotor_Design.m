function [d_welle_a, d_welle_i]=Rotor_Design(M_n, D_i)


% d_welle_a=0.35*D_a;
d_welle_a=D_i;% Shaft outer diameter corresponds to rotor inner diameter; in m

%% Design of the rotor shaft

% Material fatigue strength of the standard dimension
r_tau=0.58; % [-]
f_W=0.4; % [-]
R_mn=900; % N/mm^2
tau_WS_N=r_tau*f_W*R_mn; % N/mm^2

% technological size factor
a_dm=0.5; % [-]
d_eff_N=11; % mm
d_eff=10; % mm (corresponds to the wall thickness of the hollow shaft; is set here constantly to 10 mm to avoid iterations)
K_dm=(1-0.7686*a_dm*log10(d_eff/7.5))/(1-0.7686*a_dm*log10(d_eff_N/7.5)); % [-]

% Material fatigue strength
tau_WS=tau_WS_N*K_dm; % in N/mm^2

% Component fatigue strength
tau_AK_t=tau_WS; % in N/mm^2

% permissible torsional stress
S_Dt=2; %Safety [-]
tau_t=tau_AK_t/S_Dt; % in N/mm^2

% Minimum required torsional moment of resistance
T_max=2*M_n*1000; % in Nm
W_t=(T_max/tau_t) * 1e-6; % in m^3

% Inner diameter
% if (((d_welle_a*1000)^4)-(16*d_welle_a*1000*W_t/pi))>0
%     d_welle_i=0.001*nthroot((((d_welle_a*1000)^4)-(16*d_welle_a*1000*W_t/pi)), 4); % in meter
% else
%     d_welle_i=0;
%     disp(['Outer diameter of the shaft is only', num2str(d_welle_a*1000), ' mm'])
% end
% 
% if d_welle_a-0.02<d_welle_i
%     d_welle_i=d_welle_a-0.02;
% end
if (((d_welle_a)^4)-(16*d_welle_a*W_t/pi))>0
    d_welle_i=0.001*nthroot((((d_welle_a)^4)-(16*d_welle_a*W_t/pi)), 4); % in m
else
    d_welle_i=0;
    disp(['Outer diameter of the shaft is only', num2str(d_welle_a), ' m'])
end

if d_welle_a-0.02<d_welle_i
    d_welle_i=d_welle_a-0.02; % in m
end


