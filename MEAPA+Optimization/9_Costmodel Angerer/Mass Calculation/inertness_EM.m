function [ J ] = inertness_EM( M_nenn, n_nenn, typ_EM, i )
%in the case close to the wheel, a vector [J in the first gear J in the second gear] is created.
%regression
J_PSM = 0.0002*M_nenn-0.0029;   %[kg m^2] Pesce S.49

m_PSM = Mass_EM(M_nenn,n_nenn,'PSM');%[kg]
m_ASM = Mass_EM(M_nenn,n_nenn,'ASM');

J_ASM = (1+(m_ASM - m_PSM)/(2*m_PSM))*J_PSM;%[kg m^2]

%Lower bounds, Pesce S.49
if J_PSM < 0.001
    J_PSM = 0.001;
end
if J_ASM < 0.001
    J_ASM = 0.001;
end

tf = strcmp(typ_EM, 'PSM');
if tf == 1
    J = J_PSM*i.^2;
else
    J = J_ASM*i.^2;
end

end

