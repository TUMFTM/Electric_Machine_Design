function m_EM = Masse_SteuerEM(M_nenn, n_nenn)

P_nenn = M_nenn*n_nenn*(2*pi/60);%[W]
LD = 1800; %[W/kg] Gravimetric nominal power density from Matz
m_EM = P_nenn/LD;

end