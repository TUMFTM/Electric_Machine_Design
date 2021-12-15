% Anbauteilkosten für ein Zweistufiges Eingang Getriebe mit eTV.

%Getriebe
Oel=6;
Dichtungen=5;
Schrauben=2;
Lager=30;
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

clearvars Schaltaktuator Lageregelung Oel Dichtungen Schrauben Lager Uebriges ...
    Steuermotor_eTV Lageregelung_eTV Oel_eTV Dichtungen_eTV Schrauben_eTV Lager_eTV Uebriges_eTV