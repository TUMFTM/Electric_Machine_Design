% in this script the material cost of the engine is calculated
function [MK_gesamt, MK]=MaterialCostCalculation_EM(MV, A_wS, l_WK, typ_EM, D_a, l_fe, D, A_wL, P_n, Stkzahl, Magnetpreise, V_Gehaeuse, A_nuten_Stator, A_nuten_Rotor, A_Luftspalt, A_Magnete, D_i, Maschinendaten) %P_nenn_EM in kW
%% Materialdaten

Di_EB=7600; %density electrical sheet in kg/m^3
Di_KD=8920; %density copper wire in kg/m^3
Di_AL=2700; %density aluminum in kg/m^3
Di_ST=7700; %density steel alloy in kg/m^3
Di_MG=7400; %density magnetic material in kg/m^3

Ko_EB=1; %cost electrical sheet in €/kg
Ko_KD=10; %Cost copper wire in €/kg
Ko_AL=3.6; %Cost aluminum in €/kg
Ko_ST=1.5; %Cost steel alloy in €/kg
Ko_KU=8; %Cost of pure copper in €/kg



if Magnetpreise==1                         
    Ko_MG=81;   %Cost magnets at stable price level €/kg 
elseif Magnetpreise==2
    Ko_MG=40;   %Cost magnets at decreasing price level €/kg 
elseif Magnetpreise==3
    Ko_MG=210;  %Cost magnets at rising price level €/kg 
end


Ko_KL=0.2*P_n; %Cost in € Adhesive scaled with rated power            
Ko_HZ=0.15*P_n; %cost in € Resin scaled with rated power               
Ko_IS=0.015*P_n; %Cost in € Insulating material scaled with nominal power    
%% Material losses

Ve_EB=0.4; %punch losses in %
Ve_KD=0.1; %copper losses in %
Ve_AL=0.15; %Aluminum losses due to machining in %
Ve_ST=0.1; %steel losses due to machining in %

Stapelfaktor=0.96;

%% Cost of electrical sheet in €

                   
if strcmp(typ_EM,'ASM') && strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'Kaefiglaeufer')
v_elektroblech= (((D_a^2)/4*pi)... 
                    -A_nuten_Stator * 1e-6... 
                    -A_Luftspalt... 
                    -((D_i^2)/4*pi)... 
                    -(A_wL/1000000))...     
                    *l_fe*Stapelfaktor; %*0.7
                
elseif strcmp(typ_EM,'ASM') && strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'Schleifringlaeufer')
v_elektroblech= (((D_a^2)/4*pi)... 
                    -A_nuten_Stator * 1e-6... 
                    -A_nuten_Rotor * 1e-6...  
                    -A_Luftspalt... 
                    -((D_i^2)/4*pi))... 
                    *l_fe*Stapelfaktor; %*0.7
                
elseif strcmp(typ_EM,'PMSM') && strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'SPMSM')
v_elektroblech= (((D_a^2)/4*pi)... 
                    -A_nuten_Stator * 1e-6... 
                    -A_Luftspalt... 
                    -((D_i^2)/4*pi))... 
                    *l_fe*Stapelfaktor; %*0.7

elseif strcmp(typ_EM,'PMSM') && ((strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'IPMSM (tangential)')) || (strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'IPMSM (V-Form)')) || (strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'IPMSM (eingelassen)')))
v_elektroblech= (((D_a^2)/4*pi)... 
                    -A_nuten_Stator * 1e-6... 
                    -A_Luftspalt... 
                    -((D_i^2)/4*pi)... 
                    -A_Magnete)... 
                    *l_fe*Stapelfaktor; %r*0.7
else
    printf('Error in volume calculation for typ_EM:', typ_EM, ' and Maschienenausfuehrung:', Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung);  
end                
                
    MK_EB = v_elektroblech*Di_EB*Ko_EB*(1+Ve_EB);

%% Cost copper wire in €


if strcmp(typ_EM,'PMSM')
     v_kupferdraht = (A_wS/1000000)*(2*l_WK+l_fe); %in m^3
     MK_KD = v_kupferdraht*Di_KD*Ko_KD*(1+Ve_KD);   %in Euro
elseif strcmp(typ_EM,'ASM') && strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'Kaefiglaeufer')  
    %v_kupferdraht = ((A_wS/1000000)*(2*l_WK+l_fe)+(A_wL/1000000)*(2*l_WK+l_fe)); %in m^3
    v_kupferdraht = (A_wS/1000000)*(2*l_WK+l_fe); %in m^3
    MK_KD = v_kupferdraht*Di_KD*Ko_KD*(1+Ve_KD); %in Euro
elseif strcmp(typ_EM,'ASM') && strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'Schleifringlaeufer')
    v_kupferdraht = ((A_wS/1000000)*(2*l_WK+l_fe)+(A_wL/1000000)*(2*l_WK+l_fe)); %in m^3
    MK_KD = v_kupferdraht*Di_KD*Ko_KD*(1+Ve_KD); %in Euro
end     
    
%% Cost aluminum alloy in €
    MK_AL=V_Gehaeuse*Di_AL*Ko_AL*(1+Ve_AL); %in Euro
    
%% Cost steel alloy in €
    MK_ST=((D_a*0.55)^2*pi/4)*l_fe*2*Di_ST*Ko_ST*(1+Ve_ST);

%% Cost magnets in € (PSM only)

if strcmp(typ_EM,'PMSM')
    v_magnet = (A_Magnete)*l_fe; %Air gap min. 0.2 mm --> 0.5 mm was selected. The magnet is approx. 5 - 10 times the air gap width thick --> average value 7.5 results in a magnet thickness of 3.75 mm in kg MK_MG = v_magnet*Di_MG*Ko_MG;
% Factor 0.9: In reality, the rotor is not assembled over the complete surface 
% because certain setting points are necessary for robot gripper arms or 
% handling devices are necessary. This factor was used to 
% was approximately taken into account. In logic, the 10 % corresponds to the 
% material loss of the purchased masses, since my assumption was, 
% that material is always lost. --> but out here, because correct values
% from tool and would have to be higher instead of lower anyway.
    MK_MG =v_magnet*Di_MG*Ko_MG; % in Euro
else  
    MK_MG=0; 
end

%% Cost aluminum rotor cage in € (ASM only)

if strcmp(typ_EM,'ASM') && strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'Kaefiglaeufer')
    v_rotorkaefig=((A_wL*l_fe/1000000)+(2*D^2*pi/4*0.01)); %in m^3
    MK_RK=v_rotorkaefig*Di_AL*Ko_AL*(1+Ve_AL); %in Euro
else
    MK_RK=0; 
end

%% Cost insulation, resin and adhesive in € 

if strcmp(typ_EM,'PMSM') && (strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'SPMSM') || strcmp(Maschinendaten.Entwurf.Optionen.Maschinenausfuehrung,'IPMSM (eingelassen)'))
    MK_KL=Ko_KL+Ko_HZ+Ko_IS;   
else  
    MK_KL=Ko_HZ+Ko_IS; 
end

%% total costs

MK_gesamt=(MK_AL+MK_EB+MK_KD+MK_MG+MK_RK+MK_ST+MK_KL); % in Euro
MK.GesamtMaterial = MK_gesamt;
MK.Alulegierung = MK_AL;
MK.Elektroblech = MK_EB;
MK.Kupferdraht = MK_KD;
MK.Magnet = MK_MG;
MK.Rotorkaefig = MK_RK;
MK.Stahllegierung = MK_ST;
MK.KleberIso = MK_KL;

%% Total material costs in € incl. MGK surcharge (depending on number of pieces)

if strcmp(typ_EM,'PMSM')
    if Stkzahl>=20000
        MK_gesamt=(MK_AL+MK_EB+MK_KD+MK_MG+MK_RK+MK_ST+MK_KL)*1.1;
    elseif Stkzahl<20000 && Stkzahl>=3000
        MK_gesamt=(MK_AL+MK_EB+MK_KD+MK_MG+MK_RK+MK_ST+MK_KL)*1.16;  
    else Stkzahl<3000
        MK_gesamt=(MK_AL+MK_EB+MK_KD+MK_MG+MK_RK+MK_ST+MK_KL)*1.17;
    end
else
     if Stkzahl>=20000
        MK_gesamt=(MK_AL+MK_EB+MK_KD+MK_MG+MK_RK+MK_ST+MK_KL)*1.12;
    elseif Stkzahl<20000 && Stkzahl>=3000
        MK_gesamt=(MK_AL+MK_EB+MK_KD+MK_MG+MK_RK+MK_ST+MK_KL)*1.17;  
    else Stkzahl<3000
        MK_gesamt=(MK_AL+MK_EB+MK_KD+MK_MG+MK_RK+MK_ST+MK_KL)*1.18;
    end
end
%%
clearvars A_wL A_wS D D_a Di_AL Di_EB Di_KD Di_MG Di_ST Ko_AL Ko_EB Ko_HZ...
    Ko_IS Ko_KD Ko_KL Ko_KU Ko_MG Ko_ST l_fe MK_AL MK_EB MK_KD MK_KL MK_MG...
    MK_RK MK_ST N_S Ve_AL Ve_EB Ve_KD Ve_ST P_n
