function [Kennfeld_Darstellung,x_contourf,y_contourf]= Erstellung_Kennfeld_Darstellung (Motor)

Kennfeld_Darstellung_Antrieb=Motor.Kennfeld_Antrieb;
%Kennfeld_Darstellung_Antrieb(Kennfeld_Darstellung_Antrieb==100)=NaN;
Kennfeld_Darstellung_Antrieb=1./Kennfeld_Darstellung_Antrieb;

Kennfeld_Darstellung_Reku=Motor.Kennfeld_Reku;
%Kennfeld_Darstellung_Reku(Kennfeld_Darstellung_Reku==0)=NaN;


%Kennfeld_Darstellung=[flipud(Kennfeld_Darstellung_Reku);zeros(1,length(Kennfeld_Darstellung_Reku(1,:)));Kennfeld_Darstellung_Antrieb]

Kennfeld_Darstellung=[flipud(Kennfeld_Darstellung_Reku(2:end,:));Kennfeld_Darstellung_Antrieb]; % Annahme jedes Kennfeld hat M=0 


x_contourf=linspace(0,Motor.Drehzahlgrenze,length(Kennfeld_Darstellung(1,:)));
y_contourf=linspace(Motor.M_Max_Reku,Motor.M_Max,length(Kennfeld_Darstellung(:,1)));

end