% Anbauteilkosten f�r ein Einstufiges Eingang Getriebe mit Torque-Splitter.

%Getriebe
Oel=5;
Dichtungen=4;
Schrauben=2;
Lager=20;
Uebriges=10;

%TS
Oel_TS=6;
Dichtungen_TS=4;
Schrauben_TS=3;
Lager_TS=15;
Lamellen_TS=40;
Hydraulik_TS=20;
Uebriges_TS=10;

Anbauteile=Oel+Dichtungen+Schrauben+Lager+Uebriges+Lamellen_TS+Hydraulik_TS+Oel_TS+Dichtungen_TS+Schrauben_TS+Lager_TS+Uebriges_TS;

clearvars Oel Dichtungen Schrauben Lager Uebriges ...
Dichtungen_TS Hydraulik_TS Lager_TS Lamellen_TS Oel_TS Schrauben_TS Uebriges_TS


