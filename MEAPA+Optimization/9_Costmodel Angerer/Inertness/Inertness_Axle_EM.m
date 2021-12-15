function [ Jx, Jy, Jz ] = Inertness_Axle_EM( typ_EM, M_EM_nenn, n_EM_nenn )
%Inertias of the EM around both main axes, coefficients from
%Regression_InertiaAbout_Major_axes_EM

p_ASM_Jxz = [7.12413613020974e-05 -0.00152555980563705 0.237307431459897]; %Inertia matrix [kg * m^2]
p_ASM_Jy =  [0.000126407875744047 -0.00303300145089280 0.205162145461308]; %[kg m^2]
p_PSM_Jxz=  [5.34444530475531e-07 0.0305476020850776]; %[kg m^2]
p_PSM_Jy=   [1.18156501496769e-05 0.0257572728925912]; %[kg m^2]

if (strcmp(typ_EM,'PSM') == 1)
    Jy = p_PSM_Jy(1)*M_EM_nenn*n_EM_nenn^(2/3) + p_PSM_Jy(2); %[kg m^2]
    Jx = p_PSM_Jxz(1)*M_EM_nenn^(5/6)*n_EM_nenn^(7/6) + p_PSM_Jxz(2);
    Jz = Jx;
else
    Jy = p_ASM_Jy(1)*M_EM_nenn^2 + p_ASM_Jy(2)*M_EM_nenn + p_ASM_Jy(3);
    Jx = p_ASM_Jxz(1)*M_EM_nenn^2 + p_ASM_Jxz(2)*M_EM_nenn + p_ASM_Jxz(3);
    Jz = Jx;
end

end

 