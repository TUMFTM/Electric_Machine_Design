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

function [Maschinendaten] = Auslegung_PMSM(Maschinendaten, handles)
%   In the function Design_PMSM, you can use the literature specified to
%   design a PSMS according to the rated values. Simplifications
%   or approximations are labeled as follows
%   V/A = simplification / approximation
    
%% Variable re-storage for easier use
% Rated values
prim = Maschinendaten.Bemessungsgroessen.Primaerparameter;
sek = Maschinendaten.Bemessungsgroessen.Sekundaerparameter;
% approx. values
richt = Maschinendaten.Richtwerte;
% Options
opt = Maschinendaten.Optionen;

%% General constants
% Magnetic field constant mu_0 [H/m]
const.mu_0 = pi * 4e-7;

%%
% #########################################################################
% #   A) SIZE ESTIMATION OF ACTIVE PARTS                                  #
% #########################################################################

%%
% #########################################################################
% #   A.1) ESTIMATION OF THE MAIN DIMENSIONS                              #
% #########################################################################

%% Calculation of nominal frequency f_N [Hz] or pole pair p [-]
if(isfield(prim,'f_N'))
    prim.p = (prim.f_N * 60) / prim.n_N;
    prim.p = ceil(prim.p);
    Maschinendaten.Bemessungsgroessen.Primaerparameter.p = prim.p;
elseif(isfield(prim,'p'))
    prim.f_N = (prim.p * prim.n_N) / 60;
    Maschinendaten.Bemessungsgroessen.Primaerparameter.f_N = prim.f_N;
else
    error('Invalid input for variable frequency or number of pole pairs');
end

%% Conversion of rated power [W]
prim.P_N = prim.P_N * 1000;

%% Calculation of the synchronous rot. speed n_syn [1/min]
n_syn = (prim.f_N / prim.p) * 60;

%% Calculation of the nominal torque M_n [Nm]
M_N = prim.P_N * 60 / (2 * pi * n_syn);

%% Calculation of the phase voltage U_Str [V]
% Literature: [Mueller14, p.35 - Formula 0.6.6]
if(strcmp(sek.Schaltung,'Dreieck'))
    U_Str = prim.U_N;
elseif(strcmp(sek.Schaltung,'Stern'))
    U_Str = prim.U_N / sqrt(3);
else
    error('Invalid input');
end

%% Estimate of nominal efficiency eta_N [%] based on experience and experimental evaluation
% Literature: [Mueller08, p.571 - Fig. 9.1.7] -> Curve Fitting ->  eta_N = f(P_s)
% V/A: P_N only slightly lower than P_s -> eta_N = f(P_N)

f_eta_N = @(x) (-1.459e-07.*x.^4 + 0.0003945.*x.^3 + 98.49.*x.^2 + 1784.*x + 430.9)/(x.^2 + 18.49.*x + 5.043);
eta_N = f_eta_N(prim.P_N/1000);

%% Calculation of apparent power P_s [VA]
% Literature: [Mueller08, p.572 - Formula 9.1.17a]
P_s = prim.P_N / ((eta_N/100) * prim.cos_phi_N);

%% Calculation of phase current I_Str [A]
% Literature: [Mueller08, p.572 - Formula 9.1.17a]
I_Str = P_s / (prim.m * U_Str);

%% Calculation of nominal current I_N [A]
% Literature: [Mueller14, p.36 - Formula 0.6.7]
if(strcmp(sek.Schaltung,'Dreieck'))
    I_N = I_Str * sqrt(3);
elseif(strcmp(sek.Schaltung,'Stern'))
    I_N = I_Str;
else
    error('Invalid input');
end
 
%% Estimation of the induced voltage E_h [V] based on experience and experimental evaluation
% Literature: [Mueller08, p.603]
% V/A: Influence low, sufficient for rough design -> E_h ~ U_Str
E_h = U_Str;

% Without simplification -> E_h = f(cos_phi_N, Erregung)
    % if(prim.cos_phi_N == 1)
    %     E_h = U_Str;
    % elseif(prim.cos_phi_N == 0.8)
    %     if(strcmp(Erregung,'uebererregt'))
    %         E_h = 1.05 * U_Str;
    %     elseif(strcmp(Erregung,'untererregt'))
    %         E_h = 0.95 * U_Str;
    %     else
    %         error('Invalid input for variable excitation');
    %     end
    % else
    %     
    % end

%% Calculation of the internal apparent power P_si [VA]
% Literature: [Mueller08, p.572 - Formula 9.1.17a]
P_si = P_s * (E_h / U_Str);

%% Estimation of the utilization factor C_s [kVAmin/m^3] based on experience and experimental evaluation
% Literature: [Mueller08, p.573 - Fig. 9.1.8] -> Curve Fitting -> C_s = f(P_s/2p)

f_C_s = @(x) 5.559*exp(0.0002019*x) - 1.945*exp(-0.01073*x);
C_s = f_C_s(P_s/1000/(2*prim.p));

%% Correction of C_s depending on cos_phi_N and excitation
% Literature: [Mueller08, p.572-573]
% V/A: excitation neglected, linear correlation between given
% values for cos_phi_N (linear equation)
C_s = C_s * (0.8 + 0.25 * prim.cos_phi_N);

% Without simplification -> C_S = f(cos_phi_N, Erregung)
    % if(prim.cos_phi_N == 1)
    %     C_s = C_s * 1.05;
    % elseif(prim.cos_phi_N == 0.8)
    %     if(strcmp(Erregung,'uebererregt'))
    %         C_s = C_s;
    %     elseif(strcmp(Erregung,'untererregt'))
    %         
    %     else
    %         error('Invalid input');
    %     end
    % elseif(prim.cos_phi_N == 0)
    %     if(strcmp(Erregung,'uebererregt'))
    %         C_s = C_s * 0.8;
    %     elseif(strcmp(Erregung,'untererregt'))
    %        
    %     else
    %         error('Invalid input');
    %     end
    % else
    %   
    % end

%% Calculation of the utilisation factor C [kVAmin/m^3]
% Literature: [Mueller08, p.572 - Formula 9.1.17d]
C = C_s * (E_h / U_Str);

%% Calculation of the bore volume V_Bore [m^3]
% Literature: [Mueller08, p.572 - Formula 9.1.17c]
V_Bohrung = (P_s * 1000) / (n_syn * C_s) * (pi/4);

%% Reference value for the relative rotor length lambda [-]
% Literature: [Mueller08, p.577 - Table 9.1.3]
% if(p==1)
%     richt.lambda = 2.5; % between 1 and 4
% elseif(p>1)
%     richt.lambda = 2.0; % between 0.5 and 2.5
% end

%% Calculation of the bore diameter D_i [m]
% Literature: [Mueller08, p.576 - Formula 9.1.19a]
D_i = nthroot((P_s / 1000 * 2*prim.p) / (n_syn * C_s * richt.lambda * pi),3);
D_i_alt = nthroot((P_si / 1000 * 2*prim.p) / (n_syn * C * richt.lambda * pi),3);

% Control of the two diameters (can be commented out later)
if(D_i~=D_i_alt)
    error('Error in diameters')
end

%% Alternative calculation of bore diameter D [m]
% Literature: [Mueller08, p.576 - Formula 9.1.19b]
% Note: was not implemented

%% Check maximum circumferential speed v_max [m/s]
% Literature: [Mueller08, p.587f]
% according to [Binder18, p.209] speeds up to approx. 600km/h (~166m/s) are permissible
% according to [Meyer09, p.39] speeds up to approx. 150m/s are permissible
% V/A: Safety factor centrifugal/limit speed k_s assumed to be 3 -> maximum speed is approx. three times nominal speed

k_s = 3;
v = D_i * pi * k_s * n_syn/60;
if(v>150)
    warning('maximum circumferential speed v_max too high, adjust lambda');
end

%% Calculation of pole pitch tau_p [m]
% Literature: [Mueller08, p.577], [Mueller08, p.19 - Formula 1.1.16]
tau_p = (D_i * pi) / (2 * prim.p);

%% Calculation of the ideal rotor length l_i [m]
% Literature: [Mueller08, p.573 - Formula 9.1.17c], [Mueller08, p.577]
l_i = (P_s / 1000) / (C_s * D_i^2 * n_syn);
l_i_alt = tau_p * richt.lambda;

% Control of the two lengths (can be commented out later)
if((l_i-l_i_alt)>0.1)
    error('Error in length determination')
end

%%
% #########################################################################
% #   A.2) STATOR OUTER DIAMETER                                          #
% #########################################################################

%% Approx. value for air-gap induction B_delta [T]
% Literature: [Meyer09, p.51]
% alternative source: [Mueller08, p.582 - Table 9.1.5], [Pyr14, p.299 - Table 6.4], [Pyr14, p.299 - Table 6.4]
% if(strcmp(sek.construction_rotor,'legpole'))
% B_delta = 0.8; %between 0.8 and 1.05
% elseif(strcmp(sek.design_rotor,'full pole'))
% B_delta = 0.8; %between 0.75 and 1.05
% else
% error('Invalid input for variable construction_rotor');
% end
% 
% B_delta = 0.85; % between 0.8 and 1.05 for leg pole runners

%% Approx. value for the max. permissible induction in the yoke B_r_max [T]
% Literature: [Mueller08, p.582 - Table 9.1.5]
% if(strcmp(sek.design_rotor,'leg pole')))
% richt.B_r_max = 1.45; % between 1.0 and 1.45
% elseif(strcmp(sek.design_rotor,'full pole'))
% richt.B_r_max = 1.5; % between 1.1 and 1.5
% else
% error('Invalid input for variable construction_rotor');
% end

%% Approx. value for the max. permissible induction in the teeth B_z_max [T].
% Literature: [Mueller08, p.582 - Table 9.1.5]
% richt.B_z_max = 2.0; % between 1.6 and 2.0

%% Approx. value for the groove filling factor phi_n [-].
% Literature: [Mueller08, p.586 - Table 9.1.6]
% richt.phi_n = 0.5; % between 0.3 and 0.5 for low voltage and round wire windings

%% Approx. value for the winding factor xi_p [-].
% richt.xi_p = 0.96;

%% Calculation of the current A [A/mm]
% Literature: [Mueller08, p.572 - Formula 9.1.17d]
A = (C * sqrt(2)) / (pi^2 * richt.B_delta * richt.xi_p) * 60;

% Check whether A is in the permissible range
% Literature: [Mueller08, p.580 - Table 9.1.4]
if(strcmp(sek.Kuehlung_Stator,'Luft'))
    if(A<30 || A>120) % between 30 and 120 for indirect cooling with air
        warning('A out of bounds')
    end
elseif(strcmp(sek.Kuehlung_Stator,'Wasser'))
    if(A<160 || A>300) % between 160 and 300 for direct cooling with water
        warning('A out of bounds')
    end
else
    error('Invalid input for variable cooling');
end
% alternative Literature: [Pyr14, p.299 - Table 6.4]

% Check whether product (current coating * current density) [(A/mm)*(A/mm^2)] in permissible range
% Literature: [Mueller08, p.581 - Formula 9.1.28]
if(A*richt.S<100 || A*richt.S>350) % between 100 and 350 depending on machine size and cooling intensity
    warning('AS out of bounds')
end

%% Reference value for the current density S [A/mm^2]
% Literature: [Mueller08, p.580 - Table 9.1.4]
% if(strcmp(sec.cooling_stator,'air'))
% directional S = 3; % between 3 and 7 with indirect cooling with air 
% elseif(strcmp(sec.cooling_stator,'water'))
% directional S = 7; % between 7 and 18 with direct cooling with water 
% else
% error('Invalid input for variable cooling');
% end

% Literature: [Meyer18, p.92 - Table 5.1]
% direction S = 5; % between 3 and 6.5 for salient pole machine
%% Calculation of the minimum yoke height required h_r_min [m]
% Literature: [Meyer09, p.40 - Formula 5.13]
h_r_min = ((2/pi) * richt.B_delta * tau_p) / (2 * richt.B_r_max);

%% Calculation of the maximum yoke height h_r_max [m]
% Literature: [Meyer09, p.40 - Formula 5.14]
h_r_max = (richt.B_delta * tau_p) / (2 * richt.B_r_max);

%% Calculation of the minimum required groove height h_n_min [m]
% Literature: [Meyer09, p.41 - Formula 5.19]
h_n_min = A / (richt.S * richt.phi_n * (1- (richt.B_delta/richt.B_z_max))) / 1000;

%% Calculation of the maximum groove height required h_n_max [m]
% Literature: [Meyer09, p.41 - Formula 5.20]
h_n_max = 1.25 * h_n_min;

%% Calculation of the maximum required outer diameter D_a_max of the active material [m].
% Literature: [Meyer09, p.41 - Formula 5.21]
D_a_max = D_i + ((2.5 * A) / (richt.S * richt.phi_n * (1 - (richt.B_delta/richt.B_z_max)))) / 1000 + ((richt.B_delta * tau_p) / richt.B_r_max);

%%
% #########################################################################
% #   A.3) GEOMETRIC AIR-GAP AND ROTOR OUTER DIAMETER                     #
% #########################################################################

%% Calculation of the geometric air-gap delta_g [mm].
% Literature: [Mueller08, p.584-585 - Formula 9.1.35a-c]

% if(strcmp(sek.Bauweise_Rotor,'Schenkelpol'))
%     delta_g = 0.45 * tau_p * (A / richt.B_delta); % applies to sinusoidal field poles or arc-shaped pole contour
%     %delta_g = 0.7 * tau_p * (A / richt.B_delta); % applies to rectangular field poles
% elseif(strcmp(sek.Bauweise_Rotor,'Vollpol'))
%     delta_g = 0.25 * tau_p * (A / richt.B_delta);
% else
%     error('Ungueltige Eingabe bei Variable Bauweise_Rotor');
% end

% Literature: [Pyr14, p.305-308 - Formula 6.24a-b]

if(prim.p==1)
    delta_g = (0.2 + 0.01*prim.P_N^0.4);
else
    delta_g = (0.18 + 0.006*prim.P_N^0.4);
end
% Increase of air-gap due to variable speed operation by means of power electronics
delta_g = delta_g *1.6;

% Check if the air gap is larger than min -> reasons for production, safety
% Literature: [Mueller08, p.585 - Formula 9.1.36c]
if(delta_g<0.2)
    delta_g = 0.2;
end

%% Calculation of rotor outer diameter d_a [m]
% Literature: [Meyer09, p.42 - Formula 5.25]
d_a = D_i - 2 * (delta_g/1000);

%%
% #########################################################################
% #   A.4) AXIAL LENGTH OF ACTIVE PARTS                                   #
% #########################################################################

% V/A: Rotor and stator have the same length, partial package length
% fixed to 80mm, ventilation ducts are the same width and are only
% in the stator

%% Calculation of iron length l_e (rotor length without cooling channels) [m]
% Literature: [Mueller08, p.208-210 - Formula 2.3.25 and 2.3.30]
l_e = l_i - (2 * (delta_g/1000));

%% Approx. value for the duct width of the ventilation ducts l_v [m].
% Literature: [Meyer09, p.41], [Mueller08, p.585]
% richt.l_v = 0.01; % between 0.006 and 0.01

%% Calculation auxiliary factor gamma_v [-] (Carter factor)
% Literature: [Mueller08, p.204 - Formula 2.3.20]
% Remark: l_v in place of b_s [Mueller08, p.206 according to formula 2.3.23]
gamma_v = (1 / (1 + 5 * (delta_g/1000 / richt.l_v)));

%% Calculation of the number of required air slits n_v [-]
% Literature: [Mueller08, p.585]
if(strcmp(sek.Kuehlung_Stator,'Luft'))
    if(l_e>0.2)
        n_v = floor(l_e / 0.08); % Partial package length determined to 80mm
    else
        n_v = 0;
    end
else
    n_v = 0;
end

%% Calculation of absolute stator length l (armature length incl. cooling channels) [m]
% Literature: [Mueller08, p.208 - Formula 2.3.25]
l = l_e + (n_v * gamma_v * richt.l_v);

%%
% #########################################################################
% #   B) STATOR DESIGN                                                    #
% #########################################################################

%%
% #########################################################################
% #   B.1) Winding pattern                                                #
% #########################################################################
% Winding A1: integral slot winding, single-layer, unchorded, zoned
% Winding A2: integral slot winding, double-layer, unchorded, zoned
% Winding B: integral slot winding, double-layer, chorded, zoned
% Winding C: fractional slot winding, double-layer, chorded, zoned, Tooth coil winding (q<1) is implemented when the diameter q>1 is not allowed.

% Windingy type: predominantly chorded double-layer winding (integral or fractional slot winding) [Mueller08, p.119]

%% Choose winding type [-]
% richt.Wicklung = 'C';

%% Reference value for minimum groove pitch tau_n_min [m]
% Literature: [Meyer09, p.46], [Pyr14]-> possibly smaller (up to 0.007 m)
% richt.tau_n_min = 0.01; % between 0.01 and 0.07

%% Approx. value for the minimum number of slots per pole and phase q_min for integral slot windings [-]
% Literature: [Meyer09, p.45]
q_min = 2;

%% Calculation of the maximum stator number N_max [-]
% Literature: [Mueller08, p.19 - Formula 1.1.14], [Meyer09, p.46 - Formula 5.37]
N_max = floor((D_i * pi) / richt.tau_n_min);

%% Calculation of the minimum stator slot unmber N_min for integral slot winding [-]
% Literature: [Mueller08, p.21 - Formula 1.2.2]
N_min = (2 * prim.p * prim.m * q_min);

%% Check whether integral slot winding is possible, otherwise user input is overwritten, because a fractional slot winding must be used for realization!
if(N_min>N_max && ~strcmp(richt.Wicklung,'C'))
    warning('Winding type A and B not possible! Winding type C is used!');
    richt.Wicklung = 'C';
    Maschinendaten.Richtwerte.Wicklung = richt.Wicklung;
end

%% Determination of possíble slot numbers N_vec [-]
% Literature: self-developed, [Meyer09, p.118 - Table 4.2] for m=3

% Vector for number of grooves N_vec --> 1st boundary condition
N_vec = 1:N_max;

% Calculation of the number of slots per pole and phase vector q_vec [-]
% Literature: [Mueller08, p.21 - Formula 1.2.2], [Meyer09, p.45 - Formula 5.35]
q_vec = N_vec / (2 * prim.p * prim.m);

% Determine nominator and denominator of number of slots per pole and phase
[q_Z_vec, q_N_vec] = numden(sym(q_vec));
q_Z_vec = double(q_Z_vec);
q_N_vec = double(q_N_vec);

% Calculation of the number of slots per pole Q_vec [-]
Q_vec = N_vec / (2 * prim.p);

% Determine largest common divider of number of slots per pole and phase denominator and number of strands
for i = 1:length(q_N_vec)
    ggT_vec(1,i) = ggT_fun(q_N_vec(1,i),prim.m);
end

% Symmetry --> 2nd boundary condition
% Literature: [Mueller08, p.33 - Table 1.2.3]
% 2p/n is a natural number & ggT{n,m} = 1 & N/m is a natural number
% Note: for single layer winding the symmetry condition is tightened.
% (p/n is a natural number), but since only integral slot windings are used here as
% single layer and integral slot windings which
% automatically satisfy symmetry conditions [Mueller08, p.30]
% Further condition: Literature: [Meyer09, p.45]
% Harmonic suppression: q~=1
N_q_mat(1,:) = N_vec(~mod((2*prim.p),q_N_vec) & ~mod(N_vec,prim.m) & ggT_vec==1 & q_vec~=1);
N_q_mat(2,:) = q_vec(~mod((2*prim.p),q_N_vec) & ~mod(N_vec,prim.m) & ggT_vec==1 & q_vec~=1);
N_q_mat(3,:) = q_Z_vec(~mod((2*prim.p),q_N_vec) & ~mod(N_vec,prim.m) & ggT_vec==1 & q_vec~=1);
N_q_mat(4,:) = q_N_vec(~mod((2*prim.p),q_N_vec) & ~mod(N_vec,prim.m) & ggT_vec==1 & q_vec~=1);
N_q_mat(5,:) = Q_vec(~mod((2*prim.p),q_N_vec) & ~mod(N_vec,prim.m) & ggT_vec==1 & q_vec~=1);

% integral slot winding vs. fractional slot winding --> 3. boundary condition
switch richt.Wicklung
    case {'A1', 'A2', 'B'}
        % search for possible slot numbers N_pos for integer number of
        % slots per pole and phase
        % q_mat(p,:,3)==1 -> where denominator equals 1, is integer
        N_q_mat = N_q_mat(:,N_q_mat(4,:)==1);
    case 'C'
        % search for possible slot numbers N_pos for non-integer number of
        % slots per pole and phase
        % q_mat(p,:,3)~=1 -> where denominator is not equal to 1, is not integer
        N_q_mat = N_q_mat(:,N_q_mat(4,:)~=1);
    otherwise
        error('Invalid input for winding');
end

% feasibility of the chording--> 4th boundary condition
switch richt.Wicklung
    case {'A1', 'A2'}
        % Chording
        W_sp_rel = 1;
    case{'B', 'C'}
        % Reference value for the relative coil width W_sp_rel = W_sp / tau_p [-]
        % Literature: [Meyer09, p.47]
        % V/A: 5/6 chording for optimum harmonic suppression
        W_sp_rel = 5/6;
        
        % Sorting out the slot numbers that do not fit the chording
        if(any(~mod(N_q_mat(5,:)*W_sp_rel,1)))
            N_q_mat = N_q_mat(:,~mod(N_q_mat(5,:)*W_sp_rel,1));
        elseif(any(~mod(N_q_mat(5,:)*2/3,1)))
            % if no match was found the 5/6 chording,
            % a 2/3 chording ist attempted
            W_sp_rel = 2/3;
            N_q_mat = N_q_mat(:,~mod(N_q_mat(5,:)*W_sp_rel,1));
        else
            error('Winding type not possible!');
        end
    otherwise
        error('Invalid input for variable winding');
end
%% Selection of slot number N [-]
if(max(N_q_mat(2,:))<1) % Tooth coil winding
    % Further condition for the slot number (taking into account the winding factor of the main shaft)
    % Literature: [Mueller08, p.76]
    N_q_mat = N_q_mat(:,N_q_mat(1,:)<=(3*prim.p) & N_q_mat(1,:)>=(1.5*prim.p));

    % if a tooth coil winding (q<1) is required, then q should be
    % as close to 1/m as possible
    % Literature: [Binder12, p.94]
    [~,idx] = min(abs((1/prim.m) - N_q_mat(2,:)));
    N = N_q_mat(1,idx);
else
    % large number of slots per pole and phase required for optimum harmonic suppression -> N must be maximum
    % Literature: [Meyer09, p.45-46]
    N = max(N_q_mat(1,:));
end

% Calculation of number of slots per pole and pahse q [-]
q = N_q_mat(2,N_q_mat(1,:)==N);
q_Z = N_q_mat(3,N_q_mat(1,:)==N);
q_N = N_q_mat(4,N_q_mat(1,:)==N);

% Calculation of the groove pitch tau_n [m]
% Literature: [Meyer09, p.46]
tau_n = (D_i * pi) / N;

% Calculation of the coil width W_sp [m]
% Literature: [Meyer09, p.47]
W_sp = tau_p * W_sp_rel;

% Calculation of the groove angle alpha_n [rad]
alpha_n = (2 * pi) / N;

%% Calculation of harmonics nu [-] occurring
% Literature: [Binder12, p.124 - Formula 3.2-39]
g = [0,1:6,-1:-1:-6];
nu = 1 + ((2 .* prim.m .* g) / q_N);
nu = sort(nu,'ComparisonMethod','abs');

%% Preliminary calculations Zone factor
% Literature: [Binder12, p.128-135]
if(mod(q_N,2)) %odd
    p_u = q_N;
    Q_u = 2 * prim.m * q_Z;
else %even
    p_u = q_N / 2;
    Q_u = prim.m * q_Z;
end

if(mod(Q_u,2)) %odd
    q_1 = (Q_u + prim.m) / (2 * prim.m);
    q_2 = q_1 - 1;
else %even
    q_1 = Q_u / (2 * prim.m);
    q_2 = q_1;
end

alpha_Q = (2 * pi * p_u) / Q_u;

g_min = 0;
Y = ((g_min * Q_u) + 1) / p_u;
while(mod(Y,1))
    g_min = g_min + 1;
    Y = ((g_min * Q_u) + 1) / p_u;
end

%% Winding factor calculation
% Note: the calculation for the harmonics serves as a cross check for experts
for i = 1:length(nu)
    % Calculation of the chording factor xi_s [-]
    % Literature: [Binder12, p.136 - Formula 3.3-25]
    xi_s(i) = sin(nu(i) * (pi/2) * W_sp_rel);

    % calculation of zoning factor xi_z [-]
    % Literature: [Binder12, p.136 - Formula 3.3-25]
    xi_z(i) = (sin(nu(i) * alpha_Q * Y * (q_1/2)) - cos(nu(i) * p_u * pi * Y) * sin(nu(i) * alpha_Q * Y * (q_2/2))) / ((q_1 + q_2) * sin(nu(i) * alpha_Q * (Y/2)));

    % Calculation of the product from strain and zoning factor xi_sz [-].
    % Literature: [Meyer09, p.50]
    xi_sz(i) = xi_s(i) * xi_z(i);
end

% Winding factor for working shaft
xi_sz_haupt = abs(xi_sz(abs(nu)==1));

%%
% #########################################################################
% #   B.2) NUMBER OF COIL TURNS IN SERIES FOR A PHASE WINDING             #
% #########################################################################

%% Calculation of fundamental wave flux Phi_h [Wb]
% Literature: [Mueller08, p.199 - Formula 2.3.11]
Phi_h = (2/pi) * richt.B_delta * tau_p * l_i;

%% calculation of optimal number of coil turns in series for a phase winding w_Str_opt [-]
% Literature: [Mueller08, p.114 - Formula 1.2.89]
w_Str_opt = (sqrt(2) * E_h) / (2 * pi * prim.f_N * xi_sz_haupt * Phi_h);

%% Calculating the optimum number of parallel circuits
% Calculation of the maximum number of parallel connections a_max [-]
% and the possible combinations of parallel circuits a_pos [-] ->
% followed by selection of the best parallel connection
% best = a as small and w_Str as close to w_Str_opt as possible

% Possible combinations of parallel circuits
switch richt.Wicklung
    case 'A1'
        % max number of parallel circuits, results from condition: p/a is a natural number
        a_max = prim.p;
        
        % Possible combinations of parallel circuits
        % Condition: p/a is a natural number
        % Literature: [Mueller08, p.115]
        a_pos = 1;
        for a = 2:a_max
            if(~mod(prim.p/a,1))
                a_pos = [a_pos, a];
            end
        end
                
    case {'A2', 'B'}
        % max number of parallel circuits, results from condition: 2p/a is a natural number
        a_max = 2 * prim.p;
        
        % Possible combinations of parallel circuits
        % condition: 2p/a is a natural number
        % Literature: [Mueller08, p.115]
        a_pos = 1;
        for a = 2:a_max
            if(~mod((2*prim.p)/a,1))
                a_pos = [a_pos, a];
            end
        end
        
    case 'C'
        % Determine largest common divisor of slot number and pole pairs
        t = ggT_fun(N,prim.p);
        
        % Possible combinations of parallel circuits
        % Condition: for N/mt is even: a_max = 2t, for N/mt is odd: a_max = t
        % Literature: [Mueller08, p.115]
        if(mod(N/(prim.m*t),2))
            a_max = t;
        else
            a_max = 2 * t;
        end
        a_pos = 1:a_max;
        
    otherwise
        error('Invalid input for variable winding');
end

% Possible integer number of conductors per groove z_n
% Literature: [Mueller08, p.115], [Mueller08, p.114]
z_n_pos = floor((w_Str_opt .* 2 .* a_pos .* prim.m) ./ N);

% special conditions for double-layer windings
% Number of conductors per slot z_n must be even, so 1 is subtracted from all odd numbers of conductors.
% Literature: [Mueller08, p.115]
switch richt.Wicklung
    case 'A1'
        % nichts
    case {'A2', 'B', 'C'}
        z_n_pos(mod(z_n_pos,2)==1) = z_n_pos(mod(z_n_pos,2)==1) - 1;
    otherwise
        error('Invalid input for variable winding');
end

% possible number of coil turns in series for a phase winding
% Literature: [Mueller08, p.114]
w_pos = (N .* z_n_pos) ./ (2 .* a_pos .* prim.m);

% Deviation from optimal number
err = w_Str_opt - w_pos;

% Select minimum deviation and adopt parameters accordingly,
% if there are several minima, the solution with the smallest
% parallel connection is selected (effort and costs)
[~, idx] = min(err);
a = a_pos(1,idx);
z_n = z_n_pos(1,idx);
w_Str = w_pos(1,idx);

%% Correction of the fundamental wave flow Phi_h [Wb]
% Literature: [Mueller08, p.114 - Formula 1.2.89]
Phi_h = (sqrt(2) * E_h) / (2 * pi * prim.f_N * xi_sz_haupt * w_Str);

%% Correction of the air gap induction B_delta [T]
% Literature: [Mueller08, p.199 - Formula 2.3.11]
richt.B_delta = Phi_h / ((2/pi) * tau_p * l_i);
Maschinendaten.Richtwerte.B_delta = richt.B_delta;

%% Correction of the current A [A/mm]
% Literature: [Mueller08, p.579 - Formula 9.1.23c]
A = (2 * w_Str * prim.m * I_Str) / (pi * D_i * 1000);

% Checking whether A is in the permissible range
% Literature: [Mueller08, p.580 - Table 9.1.4]
if(strcmp(sek.Kuehlung_Stator,'Luft'))
    if(A<30 || A>120) % between 30 and 120 for indirect cooling with air
        warning('A outside the borders')
    end
elseif(strcmp(sek.Kuehlung_Stator,'Wasser'))
    if(A<160 || A>300) % between 160 and 300 for direct cooling with water
        warning('A outside the borders')
    end
else
    error('Invalid input for variable cooling');
end
% alternative Literature: [Pyr14, p.299 - Table 6.4]

% Check whether product (linear current densities * current density) [(A/mm)*(A/mm^2)] is within permissible range
% Literature: [Mueller08, p.581 - Formula 9.1.28]
if(A*richt.S<100 || A*richt.S>350) % between 100 and 350 depending on machine size and cooling intensity
    warning('AS outside the borders')
end

%%
% #########################################################################
% #   B.3) SLOT FORM                                                        #
% #########################################################################

%% Calculation of the branch current I_zw [A]
% Literature: [Mueller08, p.603]
I_zw = I_Str / a;

%% Calculation of the conductor cross-sectional area A_L [mm^2]
% Literature: [Mueller08, p.581 - Formula 9.1.28]
A_L = I_zw / richt.S;

%% Calculation of the diameter of the conductor d_L [mm]
d_L = sqrt((4 * A_L) / pi);

%% Calculation of the required groove cross-sectional area A_n [mm^2].
% Literature: [Mueller08, p.603], [Mueller08, p.586 - Formula 9.1.37]
A_n = (z_n * A_L) / richt.phi_n;

%% Calculation of the mean winding length of the windings [m]
% Literature: [Mueller08, p.586 - Formula 9.1.39]
l_m = 2 * (l + 1.3 * tau_p + (0.03 + 0.02*(U_Str/1000)));

%% Estimation of air gap flow Phi_delta [Wb]
% Literature: [Mueller08, p.604]
Phi_delta = Phi_h;

%% Calculation of the yoke flow Phi_rmax [Wb]
% Literature: [Mueller08, p.604]
Phi_rmax = Phi_delta / 2;

%% Approx. value for the iron filling factor phi_Fe [-].
% Literature: [Mueller08, p.604]
% richt.phi_Fe = 0.95;

%% Select whether or not to calculate the groove space [-]. 
% opt.Nutraumbilanz = 'JA';

%% Select whether generation of the groove should be animated [-].
% opt.Animate_Nut = 'N';

if(strcmp(opt.Nutraumbilanz,'JA'))
%% Iterative generation of the groove and tooth form in the stator
% V/A: so far only the trapezoidal shape has been executed as groove shape
% Literature: [Meyer09, p.53 - Table 5.3]

    % Approx. value for the slot width b_ns [mm].
    % Literature: [Meyer18, p.94 - Table 5.3]
    b_ns_data = [0.05 0.2 0.3 0.45; 1.5 3 3 5];
    b_ns_fun = polyfit(b_ns_data(1,:),b_ns_data(2,:),1);
    b_ns = polyval(b_ns_fun,D_i);

    % Approx. value for the slot width b_ns [mm]
    % Literature: [Meyer18, p.94 - Table 5.3]
    h_ns_data = [0.05 0.2 0.3 0.45; 0.5 1 1 2];
    h_ns_fun = polyfit(h_ns_data(1,:),h_ns_data(2,:),1);
    h_ns = polyval(h_ns_fun,D_i);

    % Specification of a minimum groove width b_no, b_n_u, b_n_m (top, bottom, middle) [mm]
    % Literature: [Meyer09, p.55 - Formula 5.49,5.50,5.51]
    if((0.2 * tau_n * 1000)<b_ns)
        b_n_u = b_ns;                                    
    else
        b_n_u = 0.2 * tau_n * 1000;
    end
    b_n_o = b_n_u;
    b_n_m = (b_n_o + b_n_u) / 2;
    
    % Determination of the thickness of the groove insulation d_iso [mm].
    % Literature: [Mueller]
    d_iso = 0.3;
    
    % Determination of slot wedge angle alpha_nk [rad]
    % Source: Assumption
    alpha_nk = pi/6;
    
    % calculation for height of slot wedge angle alpha_nk [rad]
    h_k = tan(alpha_nk) * 0.5 * (b_n_u - b_ns);
    
    % Calculation of the height to start of winding h_nk [mm]
    % V/A: Slot wedge assumed to h_k+0.5mm
    h_nk = h_ns + h_k + 0.5;
    
    % Specification of a minimum groove height h_n [mm].
    % Source: Assumption
    h_n = h_nk + 1;
    
    % Calculation of upper tooth width b_z_o [mm]
    b_z_o = (D_i * 1000 + 2 * h_n) * pi / N - b_n_o;

    % Calculation of bottom tooth width b_z_u [mm]
    b_z_u = (D_i * 1000 + 2 * h_ns + 2 * h_k) * pi / N - b_n_u;

    % Calculation of middle tooth width b_z_m [mm]
    b_z_m = (b_z_o + b_z_u) / 2;

    % Calculation of the yoke height h_r [m]
    h_r = (D_a_max - D_i)/2 - h_n/1000;
    
    % Plot of initial slot
    plot_Nut;

    % Calculation of the actual groove area A_n_tat [mm^2]
    % Source: area calculation for trapezoid
    A_n_tat = 0.5 * (2 * (x1+d_iso) + 2 * (x2+d_iso)) * (y2+d_iso - y1);

    % Upper limit iterations
    maxIter = 1e6;
    iter = 0;
    
    while(A_n_tat<A_n && iter<maxIter)
        % Calculation of the upper tooth induction B_z_o [T]
        % Literature: [Meyer09, p.55 - Formula 5.53]
        B_z_o = (richt.B_delta * tau_n * l_i) / (b_z_o / 1000 * richt.phi_Fe * l_e);

        % Calculation of the bottom tooth induction B_z_u [T]
        % Literature: [Meyer09, p.55 - Formula 5.55]
        B_z_u = (richt.B_delta * tau_n * l_i) / (b_z_u / 1000 * richt.phi_Fe * l_e);

        % Calculation of middle tooth induction B_z_m [T]
        % Source: [Meyer09, p.55 - Formula 5.54]
        B_z_m = (richt.B_delta * tau_n * l_i) / (b_z_m / 1000 * richt.phi_Fe * l_e);

        % Calculation of the yoke induction B_r [T].
        % Source: [Meyer09, p.55 - Formula 5.52]
        B_r = Phi_delta / (2 * h_r * richt.phi_Fe * l_e);
        
        % Slot shape algorithm
        % Literature: [Meyer09, p.56 - Figure 5.8]
        % Non-parallel flank grooves
%         if(B_z_o < richt.B_z_max || B_z_u < richt.B_z_max)
%             if(B_z_o < B_z_u)
%                 b_n_o = b_n_o + 0.1;
%             else
%                 b_n_u = b_n_u + 0.1;
%             end
%         elseif(B_r < B_z_o && B_r < B_z_u)
%             h_n = h_n + 0.1;
%         elseif(B_z_o < B_z_u)
%             b_n_o = b_n_o + 0.1;
%         else
%             b_n_u = b_n_u + 0.1;
%         end
        % parallel flank grooves
        if(B_z_o < richt.B_z_max && B_z_u < richt.B_z_max)
            b_n_o = b_n_o + 0.1;
            b_n_u = b_n_u + 0.1;
        elseif(B_r < B_z_o && B_r < B_z_u)
            h_n = h_n + 0.1;
        end

        % Update stator geometry
        b_n_m = (b_n_o + b_n_u) / 2;
        b_z_o = (D_i * 1000 + 2 * h_n) * pi / N - b_n_o;
        b_z_u = (D_i * 1000 + 2 * h_ns + 2 * h_k) * pi / N - b_n_u;
        b_z_m = (b_z_o + b_z_u) / 2;
        h_r = (D_a_max - D_i)/2 - h_n/1000;
        h_k = tan(alpha_nk) * 0.5 * (b_n_u - b_ns);
        h_nk = h_ns + h_k + 0.5;
        
        % Check if degenerated
        if b_z_o<=0 || b_z_u<=0 || b_z_m<=0 || h_r<=0
            error('Geometry is degenerated');
        end

        % Plot of the groove shape (required for surface calculation)
        plot_Nut;
        if(strcmp(opt.Animate_Nut,'JA'))
            pause(0.05)
        end

        % Calculation of the copper slot area A_Cu [mm^2].
        % Source: area calculation trapezoid
        A_n_tat = 0.5 * (2 * (x1+d_iso) + 2 * (x2+d_iso)) * (y2+d_iso - y1);
        
        iter = iter + 1;
    end

else
%% Estimation of the groove and tooth form in the stator

    % Estimation of the groove width b_n [mm]
    % Literature: [Mueller08, p.604]
    b_n = tau_n / 2 * 1000;

    % Estimate of the groove height h_n [mm]
    % Literature: [Mueller08, p.604]
    h_n = A_n / b_n;

    % Calculation of the tooth width b_z [mm]
    b_z = (D_i * 1000 + h_n) * pi / N - b_n;
    
    % Reference value for the slot width b_ns [mm]
    % Literature: [Meyer18, p.94 - Table 5.3]
    b_ns_data = [0.05 0.2 0.3 0.45; 1.5 3 3 5];
    b_ns_fun = polyfit(b_ns_data(1,:),b_ns_data(2,:),1);
    b_ns = polyval(b_ns_fun,D_i);
    
    % Calculation of the tooth induction B_z [T]
    % Literature: [Meyer09, p.55 - Formula 5.53]
    B_z = (richt.B_delta * tau_n * l_i) / (b_z / 1000 * richt.phi_Fe * l_e);
    
    % Checking whether B_z is in the permissible range
    if(B_z>richt.B_z_max)
        warning('B_z out of bounds')
    end
    
    % Calculation of the yoke induction B_r [T]
    % Literature: [Mueller08, p.581 - Formula 9.1.30]
    B_r = Phi_rmax / (l_e * richt.phi_Fe * h_r_max);

end

% Iteration check
if(iter>maxIter)
    warning('Check groove space balance necessary!')
end

%% Adjustment of the yoke height h_r to the maximum yoke induction [m]
% Literature: [Mueller08, p.581 - Formula 9.1.30]
if(B_r<richt.B_r_max)
    h_r = Phi_rmax / (l_e * richt.phi_Fe * richt.B_r_max);
end

%% Calculation of the relative yoke height h_r_rel [-]
% Literature: [Mueller08, p.586]
h_r_rel = h_r / tau_p;

%% Checking whether h_r_rel is in the permissible range
if(h_r_rel<0.2 || h_r_rel>0.4)
    warning('h_r_rel is out of bounds')
end

%% Calculation of the outer diameter of the stator D_a [m].
D_a = D_i + 2 * (h_r + h_n/1000);

% #########################################################################
% #   C) Rotor DESIGN                                                     #
% #########################################################################

% %% Estimation of the pole flux Phi_pk [Wb]
% % Literature: [Mueller08, p.604]
% Phi_pk = 1.2 * Phi_delta; %between 15% and 25% higher than air-gap flux
% 
% %% Calculation of the yoke flux Phi_j [Wb]
% % Literature: [Mueller08, p.604]
% Phi_j = Phi_pk / 2;
% 
% %% Guideline values for the pole inductor B_pk [T]
% % Literature: [Mueller08, p.582 - Table 9.1.5]
% if(strcmp(sek.Bauweise_Rotor,'Schenkelpol'))
%     B_pk_zul = 1.8; % zwischen 1.3 und 1.8
% elseif(strcmp(sek.Bauweise_Rotor,'Vollpol'))
%     B_pk_zul = 1.7; % zwischen 1.1 und 1.7
% else
%     error('Invalid input for variable design_rotor');
% end
% 
% %% Calculation of the cross-section of the pole core A_pk [m^2]
% % Literature: [Mueller08, p.581 - Formula 9.1.29]
% A_pk = Phi_pk / B_pk_zul;

%% Calculation of the Carter factor k_c [-]
% Literature: [Mueller08, p.203 - Formula 2.3.19], [Mueller08, p.204 - Formula 2.3.20]
gamma = 1 / (1 + 5 * delta_g / (b_ns));
k_c = tau_n / (tau_n - gamma * b_ns/1000);

%% Calculation of the ideal air gap length taking into account the groove delta_i [mm]
% Literature: [Mueller08, p. 209]
delta_i = k_c * delta_g;

%% ATTENTION
% the calculation of the inductances is normally done using a
% numerical method (FEM). Here, formulas for the analytical
% calculation are used, which is somewhat very restricted and 
% cannot guarantee accuracy
% Literature: [Ionel98, p.437]

%% Magnetic properties
% Source: [Mueller08, p.287 - Table 2.8.1]
% V/A: Neodymium iron boron, sintered, anisotropic
% mu_mr = 1.05;
% B_r_PM = 1.0;

%% Empirical factors
% Literature: [Ionel98, p.437]
% Longitudinal axes Saturation factor at idle [-]
k_so = 1.1; % between 1.1 and 1.3

if(strcmp(sek.Magnetanordnung,'IPMSM (1-schichtig, radial)'))
    % ideal Polbedeckungsfaktor [-]
    alpha_i = 0.9; % between 0.85 and 0.95 for IPM
    
    % PM Flux leakage factor [-]
    k_sigma = 0.9; % between 0.9 and 1 (low values for IPM, higher values for SPM)
    
    % Number of PM responsible for polar flux [-]
    k_b = 2; % 2 for IPM
    
    % Number of times the flow line crosses the PM [-]
    k_h = 1; % 1 for IPM
    
    % Ratio Amplitude to average value of air gap flow [-]
    k_vg = 4/pi; % 4/pi for IPM
    
    % Additional form factor for the transverse axis [-]
    k_h_q = 0; % 0 for IPM
    
    % Saturation factor for the longitudinal or transverse axis [-].
    k_sd = 1.3; % between 1.1 and 1.3
    k_sq = 1.3;
elseif(strcmp(sek.Magnetanordnung,'SPMSM (radial)'))
    % ideal Polbedeckungsfaktor [-]
    alpha_i = 2/3; % 2/3 for SPM
    
    % PM Flux leakage factor [-]
    k_sigma = 1.0; % between 0.9 and 1 (low values for IPM, higher values for SPM)
    
    % Number of PM responsible for polar flux [-]
    k_b = 1; % 1 for SPM
    
    % Number of times the flow line crosses the PM [-]
    k_h = 2; % 2 for SPM
    
    % Ratio Amplitude to average value of air gap flow [-]
    k_vg = pi/2; % pi/2 (or 1.1*alpha_i) for SPM
    
    % Additional form factor for the transverse axis [-]
    k_h_q = 1; % 1 for SPM
    
    % Saturation factor for the longitudinal or transverse axis [-]
    k_sd = 1.2; % between 1.1 and 1.3
    k_sq = 1.2; % same as k_sd for SPM
else
    error('Invalid input for variable magnet arrangement Rotor');
end

% Propotion of main inductance and total inductance [-].
k_phi = 1.25; % between 1.15 and 1.35

% Calculation of the magnet width (in the direction of magnetization) [m]
w_m = alpha_i * tau_p;

% Calculation of the magnet length (in the direction of magnetization) [m].
l_pm = (richt.B_delta * 2 * richt.mu_mr * k_so * (delta_i*1e-3)) / ...
    ((richt.B_r_PM - (richt.B_delta * (alpha_i / k_sigma) * ((pi * D_i) / (2 * prim.p * k_b * w_m)))) * k_h);
l_pm2 = (richt.mu_mr * (delta_i*1e-3)) / (((richt.B_r_PM*4*sin(w_m/(0.5*D_i)))/(richt.B_delta*pi))-1);

delta_i = k_c * delta_g + (l_pm*1e3/richt.mu_mr);

% Form factor for the longitudinal or transverse axis [-].
k_ad = (delta_i) / (delta_i + (k_h / 2) * (l_pm*1e3 / richt.mu_mr));
k_aq = (delta_i) / (delta_i + (k_h_q) * (l_pm*1e3 / richt.mu_mr));

%% Calculation of the unsaturated main inductance L_h [H]
% Literature: [Ionel98, p.437 - Formula 7]
L_h = ((prim.m*const.mu_0*(xi_sz_haupt*w_Str)^2*D_i*l_i)/(pi*prim.p^2*delta_i*1e-3));
% Literature: [Mueller08, p.528 - Formula 8.1.56], [Mueller08, p.524 - Formula 8.1.42a]
L_h2 = (prim.m/2) * (const.mu_0/(delta_i*1e-3)) * (2/pi) * tau_p * l_i * (4/pi) * ((w_Str*xi_sz_haupt)^2/(2*prim.p));

% c_1h = (6*const.mu_0)/(pi*prim.p^2)*xi_sz_haupt^2;
% r_delta = (D_i + delta_g*1e-3)/2;
% L_h1 = c_1h * w_Str^2 * l_e * (r_delta/((delta_g*1e-3)+(l_pm/richt.mu_mr)));
% L_h2 = c_1h * w_Str^2 * l_e * (r_delta/((delta_g*1e-3)+l_pm));
% k_1d = (pi/2)*alpha_i+0.5*sin(alpha_i*pi);
% k_2d = (pi/2)*(1-alpha_i)-0.5*sin(alpha_i*pi);
% L_dh = (2/pi)*(L_h1*k_1d + L_h2*k_2d);

%% Calculation of the main inductance of the longitudinal or transverse axis L_hd or L_hq [H]
% Literature: [Ionel98, p.437 - Formula 6]
L_hd = L_h * (k_ad / k_sd);
L_hq = L_h * (k_aq / k_sq);

%% Calculation of the linear current densities with formulas from Ionel (clearly different values than calculated above)
A_neu = 1e-3*(sqrt(2) * prim.p * k_sd * k_vg * alpha_i * richt.B_r_PM * (delta_i*1e-3 + (k_h / 2) * (l_pm / richt.mu_mr))) /...
     (const.mu_0 * k_phi * xi_sz_haupt * D_i * ((alpha_i / k_sigma) * ((pi * D_i) / (2 * prim.p * k_b * w_m)) + ((2 * richt.mu_mr * delta_i*1e-3 * k_so) / (k_h * l_pm))));

%% ATTENTION
% normally, for the calculation of the interlinked flux, a
% numerical method (FEM) is selected, here formulas for analytical
% calculation are used, which can be very restricted and
% cannot guarantee accuracy

%% Calculation of the interlinked PM flux [Vs]
% Literature: [Ionel98, p.437 - Formula 9]
%Psi_PM = k_phi * L_hd * I_N * sqrt(2);
%Psi_PM = xi_sz_haupt*w_Str*k_vg*alpha_i*richt.B_delta*D_i/prim.p*l_i;
Psi_PM = M_N / (3 * I_N * prim.p);

% #########################################################################
% #   D) POST CALCULATION                                                 #
% #########################################################################

%% Calculation of groove and tooth head scattering
% Approx. value for the relative scattering conductance lambda_z of the tooth head scattering [-].
% Literature: [Mueller08, p.324 - Fig. 3.7.2]
b_ns_delta_g = b_ns / delta_g;    
if(b_ns_delta_g<3)
    lambda_z_data = [1.15,1,0.8,0.62,0.5,0.39,0.3,0.26,0.19,0.13,0.1,0.07,0.03,0;...
                    0,0.125,0.3125,0.5,0.75,1,1.3125,1.5,1.75,2,2.25,2.5,2.75,3];
    lambda_z_fun = polyfit(lambda_z_data(2,:),lambda_z_data(1,:),4);
    lambda_z = polyval(lambda_z_fun,b_ns_delta_g);
else
    lambda_z_data = [0,-0.03,-0.05,-0.075,-0.1,-0.11,-0.12,-0.13,-0.14,-0.15,-0.16,-0.17,-0.18;...
                    3,4,5,6,7,8,9,10,11,12,13,14,16];
    lambda_z_fun = polyfit(lambda_z_data(2,:),lambda_z_data(1,:),4);
    lambda_z = polyval(lambda_z_fun,b_ns_delta_g);
end

if(strcmp(opt.Nutraumbilanz,'JA'))
    % Subdivision of the groove for the calculation of the groove dispersion
    % Literature: [Mueller08, p.324 - Fig. 3.7.3], [Mueller08, p.325]
    h_ue = h_nk - h_ns - h_k;
    h_l = h_n - h_nk - d_iso;
    b_n = b_n_u; %(b_n_o + b_n_u) / 2;
    b_k = (b_n + b_ns) / 2;
  
    switch richt.Wicklung
        case 'A1'
            % Calculation of the resulting relative groove-tooth head scattering conductance [-]
            % Literature: [Mueller08, p.326 - Formula 3.7.2a]
            lambda_nz = (h_l/(3 * b_n)) + (h_ue/b_n) + (h_k/b_k) + (h_ns/b_ns) + lambda_z; 

        case 'A2'
            % Approx. values for auxiliary factors k_1 and k_2 [-].
            % Literature: [Mueller08, p.330]
            k_1 = 1;
            k_2 = 1;

            % Distance between the two windings for two-layer developments d [mm]
            % Literature: [Mueller08, p.324 - Fig. 3.7.4b]
            % V/A: Distance corresponds approximately to the insulation thickness
            d = d_iso;

            % Calculation of the resulting relative groove-tooth head scattering conductance [-]
            % Literature: [Mueller08, p.330 - Formula 3.7.16]
            lambda_nz = k_1 * (h_l/(3 * b_n)) + k_2 * ((h_ue/b_n) + (h_k/b_k) + (h_ns/b_ns) + lambda_z) + (d/(4 * b_n));

        case {'B', 'C'}
            % Approx. values for auxiliary factors k_1 and k_2 [-].
            % Literature: [Mueller08, p.331 - Fig. 3.7.7]
            % V/A: m=3, single zone width
            k1_data =  [0.82,0.89,0.95,1; 2/3,0.8,0.9,1];
            k1_fun = polyfit(k1_data(2,:),k1_data(1,:),1);
            k_1 = polyval(k1_fun,W_sp_rel);

            k2_data =  [0.75,0.85,0.925,1; 2/3,0.8,0.9,1];
            k2_fun = polyfit(k2_data(2,:),k2_data(1,:),1);
            k_2 = polyval(k2_fun,W_sp_rel);

            % Distance between the two windings for two-layer windings d [mm]
            % Literature: [Mueller08, p.324 - Fig. 3.7.4b]
            % V/A: Distance corresponds approximately to the insulation thickness
            d = d_iso;

            % Calculation of the resulting relative groove-tooth head scattering conductance [-]
            % Literature: [Mueller08, p.330 - Formula 3.7.16]
            lambda_nz = k_1 * (h_l/(3 * b_n)) + k_2 * ((h_ue/b_n) + (h_k/b_k) + (h_ns/b_ns) + lambda_z) + (d/(4 * b_n)); 

        otherwise
            error('Invalid input for variable winding');
    end
else
    % Groove space balance = No -> no exact elaboration of the groove ->
    % Assumption: 90% of groove Conductor area and 10% of groove Air

    % Calculation of scattering conductance for conductor area [-]
    lambda_l = 0.9 * (h_n / (3*b_ns));
    
    % Scattering conductance Nutschlitz area [-]
    lambda_s = 0.1 * (h_n / b_ns);
    
    switch richt.Wicklung
        case {'A1', 'A2'}
            % Reference values for auxiliary factors k_1 and k_2 [-]
            % Literature: [Mueller08, p.330]
            k_1 = 1;
            k_2 = 1;

        case {'B', 'C'}
            % Reference values for auxiliary factors k_1 and k_2 [-].
            % Literature: [Mueller08, p.331 - Fig. 3.7.7]
            % V/A: m=3, single zone width
            k1_data =  [0.82,0.89,0.95,1; 2/3,0.8,0.9,1];
            k1_fun = polyfit(k1_data(2,:),k1_data(1,:),1);
            k_1 = polyval(k1_fun,W_sp_rel);

            k2_data =  [0.75,0.85,0.925,1; 2/3,0.8,0.9,1];
            k2_fun = polyfit(k2_data(2,:),k2_data(1,:),1);
            k_2 = polyval(k2_fun,W_sp_rel);

        otherwise
            error('Invalid input for variable winding');
    end

    % % Calculation of the resulting relative groove-tooth head scattering conductance [-]
    % Literature: [Mueller08, p.330 - Formula 3.7.16]
    lambda_nz = k_1 * lambda_l + k_2 * (lambda_s + lambda_z);     
end

% Calculation of the leakage inductance L_sigma_nz for slot and tooth head [H]
% Literature: [Mueller08, p.533 - Formula 8.1.76]
L_sigma_nz = 2 * const.mu_0 * l_i * (w_Str^2/prim.p) * (lambda_nz/q);

%% Calculation of the winding head dispersion
% Calculation of the conductor length l_w in the winding head [m].
% Literature: [Mueller08, p.586 - Formula 9.1.39]
l_w = 0.5 * (l_m - l);
%l_w = 1.3 * tau_p + (0.03 + 0.02 * (U_Str/1000));

switch richt.Wicklung
    case 'A1'
        % Approx. values for the relative scattering conductance lambda_ws of the winding head scattering [-].
        % Literature: [Mueller08, p.335 - Table 3.7.2]
        % V/A: m=3
        lambda_ws = 0.3;
        
        % Calculation of relative scattering conductance [-]
        % Literature: [Mueller08, p.534 - Formula 8.1.78]
        lambda_w = lambda_ws * (l_w / (l_i));
        
    case {'A2', 'B', 'C'}
        % Approx. values for the relative scattering conductance lambda_ws of the winding head scattering [-]
        % Literature: [Mueller08, p.335 - Table 3.7.2]
        % V/A: m=3, cylindrical winding
        lambda_ws = 0.3;
        
        % Calculation of the relative scattering conductance
        % Literature: [Mueller08, p.534 - Formula 8.1.79]
        % V/A: single zone width
        lambda_w = lambda_ws * (l_w / (2 * l_i));
    otherwise
        error('Invalid input for variable winding');
end

% Calculation of the assigned leakage inductance L_sigma_w for winding head [H]
% Literature: [Mueller08, p.534 - Formula 8.1.77]
L_sigma_w = 2 * const.mu_0 * l_i * (w_Str^2 / prim.p) * lambda_w;

%% Calculation of harmonic scattering
% Approx. value for the scattering coefficient of the harmonic scatter sigma_o [-].
% Literature: [Mueller08, p.340 - Fig. 3.7.12]
% V/A: m= 3
switch richt.Wicklung
    case 'A1' % W_sp_rel = 1
        sigma_o_data=[2 3 4 5 6 10000; 0.0282 0.014 0.0089 0.0065 0.0055 0.0022];      
        
    case {'A2', 'B', 'C'} % W_sp_rel = 5/6 or 2/3
        if(W_sp_rel==5/6)
            sigma_o_data=[2 3 4 5 6 10000; 0.0238 0.011 0.0065 0.0043 0.003 0.0001];
        else
            sigma_o_data=[2 3 4 5 6 10000; 0.0285 0.013 0.0088 0.0065 0.0052 0.0021];
        end
        
    otherwise
        error('Invalid input for winding');
end
sigma_o = interp1(sigma_o_data(1,:),sigma_o_data(2,:),q,'linear','extrap');

% Calculation of the assigned leakage inductance L_sigma_w for harmonic waves [H]
% Literature: [Mueller08, p.534 - Formula 8.1.79],
L_sigma_o = sigma_o * L_h;
   
%% Calculation of the slope scattering
% V/A: not required for further calculation, therefore not implemented

% Calculation of the scattering coefficient of the slope scattering sigma_schr [-]
sigma_schr = 0;

% Calculation of the assigned leakage inductance L_sigma_schr for inclination [H]
% Literature: [Mueller08, p.535 - Formula 8.1.81a]
% L_sigma_schr = sigma_schr * L_h;

%% Calculation of the leakage inductance L_sigma [H]
% Literature: [Mueller08, p.541]
L_sigma = L_sigma_nz + L_sigma_w + L_sigma_o;

%% Calculation of the synchronous inductance of the longitudinal or transverse axis L_d or L_q [H].
% Literature: [Mueller08, p.541 - Formula 8.1.105], [Mueller08, p.541 - Formula 8.1.106]
L_d = L_hd + L_sigma;
L_q = L_hq + L_sigma;
L_d = L_h + L_sigma;
L_q = L_h + L_sigma;

%% Calculation of specific conductivity of conductor materials
% V/A: only copper wire and aluminium wire implemented as conductor material

if(strcmp(sek.Material_Stator,'Kupferdraht'))
    % Reference value specific resistance rho_20 at 20°C for copper wire [mm^2/S*m]
    % Literature: [Mueller08, p.435 - Table 6.3.1]
    rho_20 = 1/58;

    % Temperature coefficient alpha for copper wire [1/K]
    % Literature: [Mueller08, p.435]
    alpha = 0.392 * 10e-3;
elseif(strcmp(sek.Material_Stator,'Aluminiumdraht'))
    % Reference value specific resistance rho_20 at 20°C for aluminium wire [mm^2/S*m]
    % Literature: [Mueller08, p.435 - Table 6.3.1]
    rho_20 = 1/37;

    % Temperature coefficient alpha for aluminium wire [1/K]
    % Literature: [Mueller08, p.435]
    alpha = 0.4 * 10e-3;
else
    error('Invalid input for variable Material_Stator');
end

% Specific resistance at theta [mm^2/S*m]
% Literature: [Mueller08, p.435 - Formula 6.3.2]
rho = rho_20 * (1 + alpha * (richt.theta - 20));    

% Specific conductivity kappa [S*m/mm^2]
% Literature: [Mueller08, p.434 - Formula 6.3.1]
kappa = 1/rho;

%% Calculation of the resistance of a winding strand R [Ohm]
% Literature: [Mueller08, p.437 - Formula 6.3.14]
R_s = (w_Str * l_m) / (a * kappa * A_L);

% #########################################################################
% #   END of design calculation                                           #
% #########################################################################

%% Re-saving of the calculated parameters in machine data struct
if(strcmp(opt.Nutraumbilanz,'JA'))
    Maschinendaten.Entwurf = struct('n_syn',n_syn,'M_N',M_N,'U_Str',U_Str,...
    'eta_N',eta_N,'P_s',P_s,'I_Str',I_Str,'I_N',I_N,'E_h',E_h,'P_si',P_si,...
    'C_s',C_s,'C',C,'V_Bohrung',V_Bohrung,'D_i',D_i,'tau_p',tau_p,'l_i',l_i,'A',A,'h_r_min',h_r_min,...
    'h_r_max',h_r_max,'h_n_min',h_n_min,'h_n_max',h_n_max,'D_a_max',D_a_max,...
    'delta_g',delta_g,'d_a',d_a,'n_v',n_v,'l',l,'l_e',l_e,'q_min',q_min,...
    'N_max',N_max,'N_min',N_min,'N_vec',N_vec,'q_vec',q_vec,'q_Z_vec',q_Z_vec,...
    'q_N_vec',q_N_vec,'Q_vec',Q_vec,'ggT_vec',ggT_vec,'N_q_mat',N_q_mat,...
    'W_sp_rel',W_sp_rel,'N',N,'q',q,'q_Z',q_Z,'q_N',q_N,'tau_n',tau_n,...
    'W_sp',W_sp,'alpha_n',alpha_n,'nu',nu,'xi_s',xi_s,'xi_z',xi_z,'xi_sz',xi_sz,...
    'xi_sz_haupt',xi_sz_haupt,'Phi_h',Phi_h,'w_Str_opt',w_Str_opt,'a_max',a_max,...
    'a_pos',a_pos,...
    'a',a,'z_n_pos',z_n_pos,'w_pos',w_pos,'err',err,'z_n',z_n,'w_Str',w_Str,...
    'I_zw',I_zw,'A_L',A_L,'d_L',d_L,'A_n',A_n,'l_m',l_m,'Phi_delta',Phi_delta,...
    'b_ns',b_ns,'h_ns',h_ns,'b_n_u',b_n_u,'b_n_o',b_n_o,'alpha_nk',alpha_nk,...
    'h_k',h_k,'d_iso',d_iso,'h_nk',h_nk,'h_n',h_n,'A_n_tat',A_n_tat,'b_z_o',b_z_o,...
    'b_z_u',b_z_u,'b_z_m',b_z_m,'h_r',h_r,'B_r',B_r,'B_z_o',B_z_o,'B_z_u',B_z_u,...
    'B_z_m',B_z_m,'b_n_m',b_n_m,'h_r_rel',h_r_rel,'D_a',D_a,'gamma',gamma,...
    'k_c',k_c,'delta_i',delta_i,'k_so',k_so,'alpha_i',alpha_i,...
    'k_sigma',k_sigma,'k_b',k_b,'k_h',k_h,'k_vg',k_vg,'k_h_q',k_h_q,'k_phi',k_phi,...
    'w_m',w_m,'l_pm',l_pm,'k_ad',k_ad,'k_aq',k_aq,'k_sd',k_sd,'k_sq',k_sq,...
    'L_h',L_h,'L_hd',L_hd,'L_hq',L_hq,'psi_PM',Psi_PM,'L_sigma_nz',L_sigma_nz,...
    'L_sigma_w',L_sigma_w,'L_sigma_o',L_sigma_o,'L_sigma',L_sigma,'L_d',L_d,...
    'L_q',L_q,'rho',rho,'kappa',kappa,'R_s',R_s);
else
    Maschinendaten.Entwurf = struct('n_syn',n_syn,'M_N',M_N,'U_Str',U_Str,...
    'eta_N',eta_N,'P_s',P_s,'I_Str',I_Str,'I_N',I_N,'E_h',E_h,'P_si',P_si,...
    'C_s',C_s,'C',C,'V_Bohrung',V_Bohrung,'D_i',D_i,'tau_p',tau_p,'l_i',l_i,'A',A,'h_r_min',h_r_min,...
    'h_r_max',h_r_max,'h_n_min',h_n_min,'h_n_max',h_n_max,'D_a_max',D_a_max,...
    'delta_g',delta_g,'d_a',d_a,'n_v',n_v,'l',l,'l_e',l_e,'q_min',q_min,...
    'N_max',N_max,'N_min',N_min,'N_vec',N_vec,'q_vec',q_vec,'q_Z_vec',q_Z_vec,...
    'q_N_vec',q_N_vec,'Q_vec',Q_vec,'ggT_vec',ggT_vec,'N_q_mat',N_q_mat,...
    'W_sp_rel',W_sp_rel,'N',N,'q',q,'q_Z',q_Z,'q_N',q_N,'tau_n',tau_n,...
    'W_sp',W_sp,'alpha_n',alpha_n,'nu',nu,'xi_s',xi_s,'xi_z',xi_z,'xi_sz',xi_sz,...
    'xi_sz_haupt',xi_sz_haupt,'Phi_h',Phi_h,'w_Str_opt',w_Str_opt,'a_max',a_max,...
    'a_pos',a_pos,...
    'a',a,'z_n_pos',z_n_pos,'w_pos',w_pos,'err',err,'z_n',z_n,'w_Str',w_Str,...
    'I_zw',I_zw,'A_L',A_L,'d_L',d_L,'A_n',A_n,'l_m',l_m,'Phi_delta',Phi_delta,...
    'b_n',b_n,'h_n',h_n,'b_z',b_z,'B_z',B_z,'Phi_rmax',Phi_rmax,...
    'h_r_rel',h_r_rel,'h_r',h_r,'D_a',D_a,'b_ns',b_ns,'gamma',gamma,'k_c',k_c,...
    'delta_i',delta_i,'k_so',k_so,'alpha_i',alpha_i,...
    'k_sigma',k_sigma,'k_b',k_b,'k_h',k_h,'k_vg',k_vg,'k_h_q',k_h_q,'k_phi',k_phi,...
    'w_m',w_m,'l_pm',l_pm,'k_ad',k_ad,'k_aq',k_aq,'k_sd',k_sd,'k_sq',k_sq,...
    'L_h',L_h,'L_hd',L_hd,'L_hq',L_hq,'psi_PM',Psi_PM,'L_sigma_nz',L_sigma_nz,...
    'L_sigma_w',L_sigma_w,'L_sigma_o',L_sigma_o,'L_sigma',L_sigma,'L_d',L_d,...
    'L_q',L_q,'rho',rho,'kappa',kappa,'R_s',R_s);
end

%% Export Excel
file_id = datestr(now,'yyyymmdd_HHMMSS');
Maschinendaten.Optionen.file_id = file_id;
folder_id = [file_id, '_data'];
Maschinendaten.Optionen.folder_id = folder_id;
mkdir('Ergebnisse',folder_id);
copyfile('Ergebnisse/Vorlage_Ergebnisse.xlsx',['Ergebnisse/',folder_id,'/Ergebnisse_',file_id,'.xlsx']);

ent_mod = rmfield(Maschinendaten.Entwurf,{'N_vec','q_vec','q_Z_vec','q_N_vec','Q_vec',...
    'ggT_vec','N_q_mat','nu','xi_s','xi_z','xi_sz','a_pos','z_n_pos','w_pos','err'});
prim_export = export_excel(Maschinendaten.Bemessungsgroessen.Primaerparameter);
sek_export = export_excel(Maschinendaten.Bemessungsgroessen.Sekundaerparameter);
richt_export = export_excel(Maschinendaten.Richtwerte);
entwurf_export = export_excel(ent_mod);

writetable(prim_export,['Ergebnisse/',folder_id,'/Ergebnisse_',file_id,'.xlsx'],'Sheet','Input Auslegung','Range','E2:E8','WriteVariableNames',0,'WriteRowNames',0)
writetable(sek_export,['Ergebnisse/',folder_id,'/Ergebnisse_',file_id,'.xlsx'],'Sheet','Input Auslegung','Range','E9:E12','WriteVariableNames',0,'WriteRowNames',0)
writetable(richt_export,['Ergebnisse/',folder_id,'/Ergebnisse_',file_id,'.xlsx'],'Sheet','Input Auslegung','Range','E13:E27','WriteVariableNames',0,'WriteRowNames',0)
writetable(entwurf_export,['Ergebnisse/',folder_id,'/Ergebnisse_',file_id,'.xlsx'],'Sheet','Output Auslegung','Range','D2:D102','WriteVariableNames',0,'WriteRowNames',0)

%% Feedback in GUI
handles.pushbutton_start.FontSize = 11;
handles.pushbutton_start.String = '<html><center>Maschinenentwurf erfolgreich abgeschlossen<br>Ergebnisse in Maschinendaten gespeichert';
var = msgbox({'Machine design successfully completed'},'Success','help','modal');
set(var,'Position',[337 673 259 80])
waitfor(var);

end