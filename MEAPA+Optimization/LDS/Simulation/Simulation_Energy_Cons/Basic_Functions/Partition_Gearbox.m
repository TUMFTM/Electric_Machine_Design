% Die Funktion übersetzt die Drehmoment-Drehzahl-Kombination aus V.
% V_Gearbox besteht aus V, wobei alle Einträge mit allen Gängen aus Gearbox übersetzt
% werden. V_Gearbox wird außerdem noch mit Informationen zum  jeweiligen Gang
% angereichert.
%
%Beispiel V=[2 3;30 30; x x]; Drehmoment_Range_1=[3; 0]; Drehmoment_Range_2=[1; -1]
%           | 5     7  |        % Moment
%  V =   | 20  20  |       % Drehzahl
%           | x     x  |        % Sonstige Informationen
%
% Gearbox.Gangabstufungen = [10 5 2];
% Gearbox.Wirkungsgrad = [0.97 0.98 0.99]; 
% 
%                           | 0.52 1.02  2.53  0.72  1.43  3.54|  % Moment
% V_Gearboxe =     | 200  100   40     200   100   40   |  % Drehzahl
%                           | x       x       x       x       x        x     |  % Sonstige Informationen 
%                           |1       2       3       1       2       3      | % Gang



function [V_Gearbox]= Partition_Gearbox(Gearbox,V)


NumofGears=length(Gearbox.gear_ratio);
NumofLines=length(V(:,1));


   
   for i=1:1:NumofLines
       x=ones(NumofGears,1)*V(i,:);
       V_R(i,:) = reshape(x,1,[]);
   end
   
   V_T_Acc=V_R(1,:);
   V_T_Reku=V_R(1,:);
   V_T_Acc(V_T_Acc<0)=0;
   V_T_Reku(V_T_Reku>0)=0;
   
   i=repmat(Gearbox.gear_ratio,1, length(V_R(1,:))/NumofGears);
   eta=repmat(Gearbox.eff,1, length(V_R(1,:))/NumofGears);
  
   %Calc Torqe
   V_T_Acc=V_T_Acc./eta./i;
   V_T_Reku=V_T_Reku.*eta./i;
   
   n=V_R(2,:).*i;
   
   
   gear=repmat(linspace(1,NumofGears,NumofGears),1, length(V));
   
   
   if NumofLines>2
   V_Gearbox=[V_T_Acc+V_T_Reku;n;V_R(3:NumofLines,:);gear];
   else
       V_Gearbox=[V_T_Acc+V_T_Reku;n;gear];
   end
   
  

end


    
    