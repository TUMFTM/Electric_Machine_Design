function     [C_m, eta, cos_phi]=Load_Factor(P_n, p, Maschinentyp)

    x5=[5 20 40 60 80 100 120 140 160 180 200];                               
    y5=[82.1 87.8 89.8 90.8 91.5 92 92.33 92.6 92.8 92.933 93];           
    n5=length(x5);                                                                            
    p5=polyfit(x5,y5,n5-1);                                                               
    x5n=5:0.1:200;                                                                            
    f5=polyval(p5,x5n);                                                                    

    % for p=1-3
    x6=[5 20 40 60 80 100 120 140 160 180 200];                                                 
    y6=[0.830 0.862 0.871 0.875 0.88 0.882 0.8822 0.8824 0.8826 0.8828 0.883];    
    n6=length(x6);                                                                                            
    p6=polyfit(x6,y6,n6-2);                                                                                 
    x6n=5:0.1:200;                                                                                              
    f6=polyval(p6,x6n);                                                                                     
    % for p=4-8
    x7=[10 20 40 60 80 100 120 140 160 180 200];                                               
    y7=[0.742 0.760 0.770 0.775 0.780 0.783 0.785 0.786 0.788 0.789 0.790];        
    n7=length(x7);                                                                                            
    p7=polyfit(x7,y7,n7-2);                                                                                 
    x7n=10:0.1:200;                                                                                            
    f7=polyval(p7,x7n);     

    
if strcmp(Maschinentyp, 'PSM')==1
    % Estimate efficiency
    eta=polyval(p5,P_n)/100;
    
    % cos_phi
    cos_phi=0.85;
    
    % Calculate utilization factor
    P_s=eta*P_n/cos_phi;
    C_m=log10(P_s/(2*p))*1.4+2.2; %Maximalwert
    C_m=log10(P_s/(2*p))*1.0+2.0; %Mittelwert
elseif strcmp(Maschinentyp, 'ASM')==1
    x=log10(P_n/(2*p));
    
    %   C_m=0.2*x^2+0.8*x+2.1; %Max. value f0r c_m 
    C_m=0.2125*x^2+0.7225*x+1.7075; %mittelwert for c_m
    
    % cos_phi
    if p>=1 && p<=3
        cos_phi=polyval(p6,P_n);
    elseif p>=4 && p<=8
        cos_phi=polyval(p7,P_n);
    else
        disp('ungueltige Polpaarzahl zur Berechnung des Leistungsfaktors')
    end
    eta=polyval(p5,P_n)/100;
    
end


