function [K_Fertigung_EM, Ko_FT] = Manufacturing_Cost(Stkzahl, Standort, typ_EM, P_n) 
%in this script the production costs of the engine are calculated
%P_n in kW

%% initialization of the production parameters depending on the selected number of pieces

if Standort==1      %Germany
    Ko_MA=25.5; %Hourly wage %Euro/h
    Ko_BB=3152; %Cost of building %Euro/m^2
elseif Standort==2  %USA
    Ko_MA=18; %Hourly wage %Euro/h
    Ko_BB=3000; %Cost of building %Euro/m^2
elseif Standort==3  %Romania
    Ko_MA=2.7; %Hourly wage %Euro/h
    Ko_BB=1492; %Cost of building %Euro/m^2
elseif Standort==4  %Czech Republic
    Ko_MA=6.7; %Hourly wage %Euro/h
    Ko_BB=1629; %Cost of building %Euro/m^2
end


if Stkzahl>=20000
     
    Abs_GB=25;          %Depreciation period for buildings in years
    Abs_MA=15;          %depreciation period for machines in years
    Tage=220;           %Working days per year
    Schichten=3;        %Shifts per day
    Stunden=8;          %hours per shift
    MA=9;               %Performing MA per shift (only MA who actually operate machines)
    OV_MA=0.4;          %Overhead for personnel in %.
    Flaeche=20000;      %Required floor space in m^2
    OEE=0.8;            %Machine OEE in %.
    Zins=0.02;          %Capital interest rate
    Energie=0.05;       %Energy cost share
    
   if strcmp(typ_EM,'PMSM')
    Investitionen_GS_PSM_AD;    %Investment costs for equipment
   else 
    Investitionen_GS_ASM_AD;    %Investment costs for equipment
   end
       
elseif Stkzahl<20000 && Stkzahl>=3000
    
    Abs_GB=25;          %Depreciation period for buildings in years
    Abs_MA=15;          %depreciation period for machines in years
    Tage=220;           %Working days per year
    Schichten=1;        %Shifts per day
    Stunden=8;          %hours per shift
    MA=9;               %Performing MA per shift (only MA who actually operate machines)
    OV_MA=0.3;          %Overhead for personnel in %.
    Flaeche=10000;      %Required floor space in m^2
    OEE=0.91;            %Machine OEE in %.
    Zins=0.02;          %Capital interest rate
    Energie=0.04;       %Energy cost share
    
   if strcmp(typ_EM,'PMSM')
    Investitionen_GS_PSM_AD;    %Investment costs for equipment
   else 
    Investitionen_GS_ASM_AD;    %Investment costs for equipment
   end
   
   
else Stkzahl<3000
    
    Abs_GB=25;          %Depreciation period for buildings in years
    Abs_MA=15;          %depreciation period for machines in years
    Tage=220;           %Working days per year
    Schichten=1;        %Shifts per day
    Stunden=8;          %hours per shift
    MA=6;               %Performing MA per shift (only MA who actually operate machines)
    OV_MA=0.2;          %Overhead for personnel in %.
    Flaeche=8000;      %Required floor space in m^2
    OEE=0.91;            %Machine OEE in %.
    Zins=0.02;          %Capital interest rate
    Energie=0.015;       %Energy cost share
    
   if strcmp(typ_EM,'PMSM')
    Investitionen_GS_PSM_AD;    %Investment costs for equipment
   else 
    Investitionen_GS_ASM_AD;    %Investment costs for equipment
   end
    
end

%% Calculate manufacturing cost structure

Ko_Personal=Schichten*Tage*Stunden*Ko_MA*MA*(1+OV_MA); %Personnel costs per year
Ko_Maschinen=Invest/Abs_MA;                            %Maschine costs per year
Ko_Instandhaltung=(1-OEE)*Ko_Maschinen;                %Maintenance costs per year
Ko_Kapital=(Flaeche*Ko_BB+Invest)*Zins;                %Capital costs per year
Ko_Abs=(Flaeche*Ko_BB)/Abs_GB;                         %Depreciation costs per year
Ko_Energie=Invest*Energie;                             %Energy costs per year

%Total production costs per year
Ko_Fertigung=Ko_Personal+Ko_Maschinen+Ko_Instandhaltung+Ko_Kapital+Ko_Abs+Ko_Energie;

Ko_FT.gesamt = Ko_Fertigung;
Ko_FT.personal = Ko_Personal;
Ko_FT.machine = Ko_Maschinen;
Ko_FT.instandhaltung = Ko_Instandhaltung;
Ko_FT.kapital = Ko_Kapital;
Ko_FT.abschreibung = Ko_Abs;
Ko_FT.energie = Ko_Energie;

%divided among the number of units
cpp=Ko_Fertigung/Stkzahl; % in cost per year per unit

%% Scaling of manufacturing costs to nominal output.
%Scaling of manufacturing costs to the nominal power of the machine. With
%the nominal power of the machine, the manufacturing costs increase
%logarithmically according to Ehrlenspiel [99] (design is based on 30 kW
%nominal power)

if P_n<=20
    cpp=cpp*0.9;
elseif P_n>20 && P_n<=35
    cpp=cpp;
elseif P_n>35 && P_n<=50
    cpp=cpp*1.1;
elseif P_n>50 && P_n<=60
    cpp=cpp*1.15;
elseif P_n>60 && P_n<=70
    cpp=cpp*1.2;
elseif P_n>70 && P_n<=85
    cpp=cpp*1.25;
elseif P_n>85 && P_n<=100
    cpp=cpp*1.33;
elseif P_n>100 && P_n<=130
    cpp=cpp*1.40;
elseif P_n>130
    cpp=cpp*1.48;    
end

K_Fertigung_EM=cpp; % in cost per year per unit
Ko_FT.Fertigungskosten_je_Stk = K_Fertigung_EM;
%%

clearvars Abs_GB Abs_MA Energie Flaeche Invest Ko_Abs Ko_BB Ko_Energie Ko_Fertigung...
    Ko_Instandhaltung Ko_Kapital Ko_MA Ko_Maschinen Ko_Personal MA OEE OV_MA Stunden...
    Tage Zins Schichten