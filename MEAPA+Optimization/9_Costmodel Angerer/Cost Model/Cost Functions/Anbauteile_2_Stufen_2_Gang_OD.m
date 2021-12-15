% Add-on part costs for a single-stage two-speed transmission with open
% differential.

%Gear
Oel=7;
Dichtungen=5;
Schrauben=2;
Lager=30;
Schaltaktuator=15;
Lageregelung=10;
Uebriges=12;

%Differential
Oel_OD=2;
Dichtungen_OD=2;
Schrauben_OD=2;
Lager_OD=10;
Uebriges_OD=5;

Anbauteile=Schaltaktuator+Lageregelung+Oel+Dichtungen+Schrauben+Lager+Uebriges+Oel_OD+Dichtungen_OD+Schrauben_OD+Lager_OD+Uebriges_OD;

clearvars Schaltaktuator Lageregelung Oel Dichtungen Schrauben Lager Uebriges ...
    Oel_OD Dichtungen_OD Schrauben_OD Lager_OD Uebriges_OD