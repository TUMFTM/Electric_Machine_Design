function [y, cons, dimensions] = optimization_Kalt_PMSM_objfun(x)
% Objective function : Test problem 'CONSTR'.
%*************************************************************************

y = [0,0,0]; %Output, angepasst auf 3 Variablen
cons = [0,0,0]; %Constrains, z.B. Strom auf bestimmte Ampere begrenzen


%Maschinenauslegung
% Maschinentyp - ASM oder PMSM
rated.type = x(5);
% Nennleistung P_N [W]
% test = 5:5:200;
% rated.P_N   = test(x(1)); --> wenn z.B. nur 5er schritte gewünscht-->upper und lower bounds anpassen!
rated.P_N   = x(1)*1000;
% Nenndrehzahl n_N [U/min]
rated.n_N   = x(2);
% Nennspannung U_N [V]
rated.U_N   = x(3);
% Maximal Drehzahl n_max [U/min]
rated.nmax = 12000;
% Polpaarzahl p [-]
rated.p     = x(4);
% Nennfrequenz f_N [-]
rated.f_N   = (rated.p * rated.n_N) / 60;
% Strangzahl m [-]
rated.m     = 3;
% Leistungsfaktor [-]
rated.cos_phi_N = 0.95;

%Kühlung variabel
rated.Kuehlungsart = x(6);

%Magnetanordnung variabel
rated.Magnetanordnung = x(7);

% Set to default BMW i3 vehicle
rated.LDS.fz_m.String = num2str(1640);
rated.LDS.cW.String = num2str(0.3);
rated.LDS.A.String = num2str(2.38);
rated.LDS.tyre_r.String = num2str(0.3498);
rated.LDS.battcap.String = num2str(x(8)); %num2str(22); %Batterie
rated.LDS.aux.String = num2str(1500);
rated.LDS.GearRatio.String = num2str(x(9)); %num2str(9.7); %Gear

%Anzahl Maschinen
rated.LDS.AnzMasch.String = num2str(x(10));

if rated.LDS.AnzMasch.String == 1
    rated.LDS.AnzMasch.String = 'GM_X';
else
    rated.LDS.AnzMasch.String = 'GM_GM';
end

rated.LDS.visual_LDS = 0; %visual of LDS desired: 1, visual not desired: 0

[Entwurf, Analyse, Optimierer] = MEAPA_Skript(rated);

y(1) = Optimierer.Kosten_in_Euro;
y(2) = - Optimierer.Zykluseffizienz_in_Prozent; % Minus, da maximiert werden soll
y(3) = Optimierer.Volumen_in_m3;

%dimensions = Entwurf.Geometrie.l_Fe; %hier beliebig anpassbare Variable
dimensions = Entwurf;

end

