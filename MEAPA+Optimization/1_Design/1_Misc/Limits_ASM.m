% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Modell for the design and analysis of PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (svenja.kalt@tum.de), 
%           Jonathan Erhard 
% -------------------------------------------------------------------------

% Rated power per pole pair in W
% pos.P_N = [max min (p==8); max min (p==6); ... ; max min (p==1)]
limit.P_N = [1600e3 64e3; 5000e3 1.0e3; 5000e3 1.0e3; 5000e3 1.0e3; 5000e3 1.0e3; 4000e3 1.0e3; 2000e3 1.2e3];

% Rated rotational velocity in 1/min
limit.n_N = [Inf 0];

% Rated voltage in V
limit.U_N = [Inf 0];

% Number of pole pairs
limit.p = [8 6 5 4 3 2 1];

% Rated Frequency in Hz
limit.f_N = [Inf 0];

% Number of bars
limit.m = [3 3];