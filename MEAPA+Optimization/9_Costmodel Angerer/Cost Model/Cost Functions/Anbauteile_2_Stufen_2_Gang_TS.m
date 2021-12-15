% Add-on part costs for a two-stage multi-speed transmission with torque splitter.

%Gear
Oel=7;
Dichtungen=5;
Schrauben=2;
Lager=30;
Schaltaktuator=15;
Lageregelung=10;
Uebriges=12;

%TS
Oel_TS=6;
Dichtungen_TS=4;
Schrauben_TS=3;
Lager_TS=15;
Lamellen_TS=40;
Hydraulik_TS=20;
Uebriges_TS=10;

Anbauteile=Schaltaktuator+Lageregelung+Oel+Dichtungen+Schrauben+Lager+Uebriges+Lamellen_TS+Hydraulik_TS+Oel_TS+Dichtungen_TS+Schrauben_TS+Lager_TS+Uebriges_TS;

clearvars Oel Dichtungen Schrauben Lager Uebriges Schaltaktuator Lageregelung ...
Dichtungen_TS Hydraulik_TS Lager_TS Lamellen_TS Oel_TS Schrauben_TS Uebriges_TS


