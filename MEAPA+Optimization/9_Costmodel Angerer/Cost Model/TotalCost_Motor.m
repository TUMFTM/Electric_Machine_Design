%model was analytically %designed for quantities between 2,000 and 100,000.
%design. Here, parametric projection according to Ehrlenspiel [99, p. 175] 
%for unit numbers up to 2,000,000 whereby a conservative classification was made.

Kosten_Motor=MK_gesamt+K_Anbauteile+K_Fertigung_EM; 

if Stkzahl<=100000
    Kosten_Motor=Kosten_Motor;
elseif Stkzahl>100000 && Stkzahl<=250000
    x=(1-0.92)/(100001-250000);
    y=1-x*(100001-Stkzahl);
    Kosten_Motor=Kosten_Motor*y;
elseif Stkzahl>250000 && Stkzahl<=500000
    Kosten_Motor=Kosten_Motor*0.92;
    x=(1-0.97)/(250001-500000);
    y=1-x*(250001-Stkzahl);
    Kosten_Motor=Kosten_Motor*y;
elseif Stkzahl>500000
    Kosten_Motor=Kosten_Motor*0.8924;
    x=(1-0.92)/(500001-2000000);
    y=1-x*(500001-Stkzahl);
    Kosten_Motor=Kosten_Motor*y;
end


RES.em.K_ges=Kosten_Motor;