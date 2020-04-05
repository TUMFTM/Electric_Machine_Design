% -------------------------------------------------------------------------
% TU Muenchen - Lehrstuhl fuer Fahrzeugtechnik (FTM)
% -------------------------------------------------------------------------
% Modell fuer den Entwurf und die Analyse einer PMSM oder ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor: Svenja Kalt (kalt@ftm.mw.tum.de)
%        Jonathan Erhard
% -------------------------------------------------------------------------

% Hinweis:

% Nennleistung je Polpaar in W
% pos.P_N = [max min (p==8); max min (p==6); ... ; max min (p==1)]
limit.P_N = [1600e3 64e3; 5000e3 1.0e3; 5000e3 1.0e3; 5000e3 1.0e3; 5000e3 1.0e3; 4000e3 1.0e3; 2000e3 1.2e3];

% Nenndrehzahl in 1/min
limit.n_N = [Inf 0];

% Nennspannung in V
limit.U_N = [Inf 0];

% Polpaarzahl
limit.p = [8 6 5 4 3 2 1];

% Nennfrequenz in Hz
limit.f_N = [Inf 0];

% Stangzahl
limit.m = [3 3];
