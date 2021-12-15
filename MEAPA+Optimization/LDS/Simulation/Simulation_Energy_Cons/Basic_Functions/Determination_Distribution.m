function [P_EL_MIN, DATA] = Determination_Distribution (P_EL_1,P_EL_2, DATA_IN)

NumofLines=length(P_EL_1(:,1));

P_EL_ALL=P_EL_1(1,:)+P_EL_2(1,:);

n=100;

P_EL_ALL_R=reshape(P_EL_ALL,n,[]);

[~, idx]=min(P_EL_ALL_R,[ ],1);

V_idx=([0:1:P_EL_1(end,end)-1]*n)+idx;

if NumofLines>2
P_EL_ALL(2:NumofLines-1,:)=P_EL_1(2:NumofLines-1,:);
end

P_EL_MIN=P_EL_ALL(:,V_idx);

DATA=DATA_IN(:,V_idx);

  
end