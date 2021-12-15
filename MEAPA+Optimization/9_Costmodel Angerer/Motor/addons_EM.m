function K_Anbauteile=addons_EM(Stkzahl, P_n, akt_az)
%In this script, the attachment costs for an electric motor are calculated.
% The selection of the mounting parts is variable with the number of units and the bearings are
% scale with the rated power (from SKF bearing catalog).
% Source for unit-dependent mounting part costs of electric motors [12].

if akt_az==1
    P_n_VA_ANB=P_n; % in kW
else 
    P_n_VA_ANB=P_n*2;
end


if Stkzahl>=20000
    Resolver=29;
    Stecker=3;
    NTC=0.8;
    Schrauben=1+0.0111*P_n_VA_ANB;
    Lager=8+0.133*P_n_VA_ANB;
    Uebriges=3+0.02*P_n_VA_ANB;
    
    K_Anbauteile=Resolver+Stecker+NTC+Schrauben+Lager+Uebriges; % in Euro per piece
      
elseif Stkzahl<20000 && Stkzahl>=3000
    Resolver=32;
    Stecker=4;
    NTC=1;
    Schrauben=3;
    Lager=20+0.333*P_n_VA_ANB;
    Uebriges=5;
    
    K_Anbauteile=Resolver+Stecker+NTC+Schrauben+Lager+Uebriges; % in Euro per piece
      
else Stkzahl<3000
    Resolver=40;
    Stecker=10;
    NTC=2;
    Schrauben=3;
    Lager=40+0.333*P_n_VA_ANB;
    Uebriges=7;
    
    K_Anbauteile=Resolver+Stecker+NTC+Schrauben+Lager+Uebriges; % in Euro per piece
end
  %%  
clearvars Resolver Stecker NTC Schrauben Lager Uebriges P_n_VA_ANB