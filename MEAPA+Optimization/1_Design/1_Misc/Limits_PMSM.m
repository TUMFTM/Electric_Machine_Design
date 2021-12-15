% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Modell for the design and analysis of PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (svenja.kalt@tum.de), 
%           Jonathan Erhard 
% -------------------------------------------------------------------------


% Rated power in kW
limit.P_N = [2000e3 0];

% Rated rotational speed in 1/min
limit.n_N = [Inf 0];

% Rated voltage in V
limit.U_N = [Inf 0];

% number of pole pairs
limit.p = [Inf 0];

% Rated frequency in Hz
limit.f_N = [Inf 0];

% Power factor
limit.cos_phi_N = [1.0 0.7];

% Number of bars
limit.m = [3 3];