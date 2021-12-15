function [P_EL_Gear, DATA] = Determination_Gear (P_EL, Gearbox, DATA_IN)

NumberofGears= length(Gearbox.gear_ratio);

NumofLines=length(P_EL(:,1));

P_EL_R=reshape(P_EL(1,:),NumberofGears,[]);

[~, idx]=min(P_EL_R,[ ],1);

V_idx=[0:NumberofGears:(length(P_EL(1,:))-NumberofGears)]+idx;

P_EL_Gear=P_EL(1:NumofLines-1,V_idx);

DATA=[DATA_IN(:,V_idx);P_EL(NumofLines,V_idx)];































% 
% NumberofGears= length(Gearbox.gear_ratio);
% 
% NumbofLines=length(P_EL(:,1));
% %P_EL_Gear=[];
% 
% P_EL_Gear=zeros(NumbofLines-1,length(P_EL(1,:))/NumberofGears);
% 
%     j=1;
%     
%   %   GearPossibilities=[];
%    
%  for i=1:NumberofGears:length(P_EL(1,:))   
%      GearPossibilities=[];   
%      for x=i:1:(i+NumberofGears-1)
%      GearPossibilities=[GearPossibilities P_EL(1,x)];
%      end
% 
%      [P_EL_G(1,j), Gear]=min(GearPossibilities);
%         
%     
%      
%      DATA(:,j)=[DATA_IN(:,i+Gear-1);Gear];
% %     if i==1
%      
%      if NumbofLines-1 >=2
%      P_EL_Gear(:,j) = [P_EL_G(:,j);P_EL(2:NumbofLines-1,x)];  
%      else
%      P_EL_Gear(:,j) = [P_EL_G(:,j)];
%      end    
%      
         
%     else

%     if NumbofLines-1 >=2
%     P_EL_Griebe = [P_EL_Griebe [P_EL_G(:,j);P_EL(2:NumbofLines-1,x)]];  
%     else
%     P_EL_Griebe = [P_EL_Griebe [P_EL_G(:,j)]];
%     end
     
%    end
     
     
     
     
%      
%      j=j+1;
%      
%  end

  %for i=1:NumberofGears:length(P_EL(1,:))   
  %   [P_EL(j) idx_Gang(j)]=min([P_EL(i) P_EL(i+1) P_EL(i+2)]);
  %      
  %   j=j+1;
         end
 
 
 

