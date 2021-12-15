
function [P_EL_Motor, DATA] = Determination_Power (V, Motor, vehicle)

NumofLines=length(V(:,1));

DATA=zeros(2,length(V(1,:)));
P_EL_Motor=zeros(NumofLines-1,length(V(1,:)));

[xd,yd]=ndgrid(Motor.eff_n_axis,Motor.eff_T_axis);

GRIDINT=griddedInterpolant(xd,yd, Motor.eff','linear','none'); %Objekt zur Interpolation, Lookup-Table
%GRIDINT=griddedInterpolant(Motor.eff_n_axis,Motor.eff_T_axis, Motor.eff','linear','none'); %Objekt zur Interpolation, Lookup-Table

%EFF=GRIDINT(V(2,:).*60, V(1,:)); %Übeltaeter
EFF=GRIDINT(V(2,:).*60, V(1,:)); 

P_EL_Motor(1,:)=V(1,:).*V(2,:).*2.*pi./EFF;

if NumofLines>2
    P_EL_Motor(2:NumofLines-1,:)=V(3:NumofLines,:);
end

DATA(:,:)=V(1:2,:); % DATA-Matrix mit Moment und Drehzahl zu den einzelnen Leistungen zur Nachvollziehung

end