% -------------------------------------------------------------------------
% TU Muenchen - Lehrstuhl fuer Fahrzeugtechnik (FTM)
% -------------------------------------------------------------------------
% Modell fuer den Entwurf und die Analyse einer PMSM oder ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Jonathan Erhard (jonathan.erhard@tum.de)
% Version:  v0.7 - 19.05.2019 - je
% -------------------------------------------------------------------------

% Hinweis:

% Nennleistung in kW
limit.P_N = [2000e3 0];

% Nenndrehzahl in 1/min
limit.n_N = [Inf 0];

% Nennspannung in V
limit.U_N = [Inf 0];

% Polpaarzahl
limit.p = [Inf 0];

% Nennfrequenz in Hz
limit.f_N = [Inf 0];

% Leistungsfaktor
limit.cos_phi_N = [1.0 0.7];

% Stangzahl
limit.m = [3 3];