% Funktion erstellt ausgehend von der Eingangsmatrix v und den erhlaubten
% Drehmomentenspreizungen Torque_Range_1 und Torque_Range_2 
% zwei zusammengehörige Matritzen V_Split_1 und V_Split_2, die die möglichen
% Drehmomentenkombinationen darstellen, um V zu erfüllen.
% V_Split_1(1,i)+V_Split_2(1,i)=V(1,?)
% 
% Beispiel V=[2;30]; Torque_Range_1=[3; 0]; Torque_Range_2=[1; -1]

% V_Split_1= |3  2  1 |    V_Split_2= |-1 0  1 |    Drehmomente
%                   |30 30 30|                    |30 30 30|    Drehzahl
%                     |1  1  1 |                     |1  1  1 |    i von V(1,i)

function [V_Split_1, V_Split_2]= Partition_Torque(Torque_Range_1, Torque_Range_2, V)

%NumofLines=length(V(:,1));


  n = 100;



V_TR11=ones(1,length(V(1,:)))*Torque_Range_1(1);
V_1_Help= V(1,:)-Torque_Range_2(2);

V_TR21=ones(1,length(V(1,:)))*Torque_Range_2(1);
V_2_Help= V(1,:)-Torque_Range_1(2);

V_TR12=ones(1,length(V(1,:)))*Torque_Range_1(2);
V_3_Help= V(1,:)-Torque_Range_2(1);

V_TR22=ones(1,length(V(1,:)))*Torque_Range_2(2);
V_4_Help= V(1,:)-Torque_Range_1(1);


Max_1=min(V_TR11,V_1_Help);
Max_2=min(V_TR21,V_2_Help);
Min_1=max(V_TR12, V_3_Help);
Min_2=max(V_TR22, V_4_Help);


Max_1(Max_1<Min_1)=0;

Max_2(Max_2<Min_2)=0;


Help=linspace(0,1,n);

V_Split_1=Min_1+Help'.*(Max_1-Min_1);

V_Split_1=reshape(V_Split_1, 1,[]);

   for i=1:1:length(V(:,1))
       x=ones(n,1)*V(i,:);
       V_R(i,:) = reshape(x,1,[]);
   end

V_Split_2=V_R(1,:)-V_Split_1;

Partition=[1:1:length(V(1,:))];
Partition=ones(n,1)*Partition;
Partition= reshape(Partition,1,[]);

V_Split_1=[V_Split_1;V_R(2:length(V(:,1)),:);Partition];

V_Split_2=[V_Split_2; V_R(2:length(V(:,1)),:);Partition];



end