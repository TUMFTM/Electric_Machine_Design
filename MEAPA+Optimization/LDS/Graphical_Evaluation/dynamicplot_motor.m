function [Tn] =dynamicplot_motor(vehicle, dc, Tn)

%% Rearange vehicle.MOTOR
 vehicle.MOTOR=vehicle.MOTOR(~cellfun('isempty',vehicle.MOTOR)) ; 

figure
%% Check amount of motors
AmountMot=length(vehicle.MOTOR);

%% Creation of maps
for i=1:1:AmountMot
    MotMap{i}=CreationMapPlot(vehicle.MOTOR{1,i}.eff);
end

%% Load Colors
 [TumColors] = tumColors(  );
 
 
 %% Framework Motors
  Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
 for i=1:1:AmountMot
subplot(round(AmountMot/2)+1,2,i)
contourf(vehicle.MOTOR{1,i}.eff_n_axis, vehicle.MOTOR{1,i}.eff_T_axis,MotMap{i}, Plot_Label_Vektor_eta) %Beim Weglassen dieser Zeile, bleibt das Kennfeld im Hintergrund weg

title(strcat('Motor',int2str(i)))
xlabel('Rev. Speed in 1/min')
ylabel('Torque in Nm')

hold on

 end

 %% Framework Cycle
subplot(round(AmountMot/2)+1,2,(2*round(AmountMot/2)+1):(2*round(AmountMot/2)+2));
plot (dc.speed)
title('Driving Cycle')
xlabel('Time in s')
ylabel('Speed in m/s')

hold on


%% First Scatter
for i=1:1:AmountMot
 subplot(round(AmountMot/2)+1,2,i)
 s{i}=scatter(Tn(2,i),Tn(1,i), 'filled');
end
subplot(round(AmountMot/2)+1,2,(2*round(AmountMot/2)+1):(2*round(AmountMot/2)+2));
s{AmountMot+1}=scatter(1, dc.speed(1),'filled');

%% Other Scatter

for i = 2:length(Tn(1,:)) 
    
    for j=1:AmountMot
    s{j}.CData=TumColors.secondaryDarkGrey;
    end
    s{AmountMot+1}.CData=TumColors.primaryBlue;
    
    
     for j=1:AmountMot
    s{j}=scatter(subplot(round(AmountMot/2)+1,2,j),Tn(j*2,i)*60,Tn(j*2-1,i), 'filled'); %Tn stellt Matrix der Betriebspunkte dar, durch *60 ist es danach in Minuten. Vorher ist es in Sekunden, jetzt Umdrehungen/Minute.
     end
       % pause 2/10 second: 
       
    s{AmountMot+1}=scatter(subplot(round(AmountMot/2)+1,2,(2*round(AmountMot/2)+1):(2*round(AmountMot/2)+2)),i, dc.speed(i),'filled');   
       
    %pause(0.0000001) %war vorher auf 0.05 %hier kann plotpause ausgestellt werden
end


end

