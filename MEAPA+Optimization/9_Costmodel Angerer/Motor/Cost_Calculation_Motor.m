function [Kosten_Motor, MK, K_Anbauteile, KO_FT]=Cost_Calculation_Motor(MV, A_wS, l_WK, typ_EM, D_a, l_fe, D, A_wL, P_n, Stkzahl, Magnetpreise, akt_az, Standort, V_Gehaeuse, A_nuten_Stator, A_nuten_Rotor, A_Luftspalt, A_Magnete, D_i, Maschinendaten)
% In this script, the individual cost items are calculated on the basis of the
% machine dimensioning is calculated.  


%% Material cost calculation

[MK_gesamt, MK]  = MaterialCostCalculation_EM(MV, A_wS, l_WK, typ_EM, D_a, l_fe, D, A_wL, P_n, Stkzahl, Magnetpreise, V_Gehaeuse, A_nuten_Stator, A_nuten_Rotor, A_Luftspalt, A_Magnete, D_i, Maschinendaten);       %% Add-on cost calculation

K_Anbauteile = addons_EM(Stkzahl, P_n, akt_az);              

%% Production cost calculation

[K_Fertigung_EM, KO_FT] = Manufacturing_Cost(Stkzahl, Standort, typ_EM, P_n);             

%% Toatl cost calculation motor

TotalCost_Motor

