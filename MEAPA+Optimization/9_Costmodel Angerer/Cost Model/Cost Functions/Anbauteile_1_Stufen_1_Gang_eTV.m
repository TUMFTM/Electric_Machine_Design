% Mounting part costs for a single-stage input transmission with eTV.

%Gear
Oel=5;
Dichtungen=4;
Schrauben=2;
Lager=20;
Uebriges=10;

%eTV
Oel_eTV=6;
Dichtungen_eTV=5;
Schrauben_eTV=4;
Lager_eTV=30;
Steuermotor_eTV=90;
Lageregelung_eTV=10;
Uebriges_eTV=35;

Anbauteile=Oel+Dichtungen+Schrauben+Lager+Uebriges+Steuermotor_eTV+Lageregelung_eTV+Oel_eTV+Dichtungen_eTV+Schrauben_eTV+Lager_eTV+Uebriges_eTV;

clearvars Oel Dichtungen Schrauben Lager Uebriges ...
    Steuermotor_eTV Lageregelung_eTV Oel_eTV Dichtungen_eTV Schrauben_eTV Lager_eTV Uebriges_eTV


