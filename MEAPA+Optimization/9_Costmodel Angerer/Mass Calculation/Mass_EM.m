function m_EM = Mass_EM(M_nenn, n_nenn, typ_EM)

% Pesce
P_nenn = M_nenn*n_nenn*(2*pi/60);%[W]
compPSM = strcmp('PSM', typ_EM);
if compPSM == 1
    m_PSM = 0.5*((13.847*log(P_nenn/1000)-13.003)+(10.979*log(M_nenn)-17.908));
    m_EM = m_PSM;
end
compASM = strcmp('ASM', typ_EM);
if compASM == 1
    m_ASM = 0.5*((16.152*log(P_nenn/1000)+3.4576)+(24.184*log(M_nenn)-49.245));
    m_EM = m_ASM;
end

end