function [] dynamicplot_motor(vehicle, dc, Tn)

%% Check amount of motors
AmountMot=length(vehicle.MOTOR)

%% Creation of maps
for i=1:1:AmountMot
    MotMap{i}=CreationMapPlot(vehicle.MOTOR{1,i}.eff)
end

%% Load Colors
 [TumColors] = tumColors(  );
 
 
 %% Framework Motors
 
 for i=1:1:AmountMot
subplot(round(AmountMot/2)+1,2,i)
contourf(vehicle.MOTOR{1,i}.eff_n_axis, vehicle.MOTOR{1,i}.eff_T_axis,MotMap{i})

title(strcat('Motor',int2str(i)))

hold on

 end

 %% Framework Cycle
subplot(round(AmountMot/2)+1,2,i+1:i+2);
plot (dc.speed)
title('Driving Cycle')


%% First Scatter
for i=1:1:AmountMot+1
 subplot(round(AmountMot/2)+1,2,i)
 s{i}=scatter(Tn(2,i),Tn(1,i), 'filled');
end

%% Other Scatter

for i = 2:length(Tn(1,:)) 
    
    for j=1:AmountMot-1
    s1.CData=TumColors.secondaryDarkGrey;
    s2.CData=TumColors.secondaryDarkGrey;
    s3.CData=TumColors.primaryBlue;
    
     s1=scatter(subplot(2,2,1),Drehzahl_1(i),Moment_1(i), 'filled');
     
     s2=scatter(subplot(2,2,2),Drehzahl_2(i),Moment_2(i), 'filled');
     
     s3=scatter(subplot(2,2,3:4),i, Zyklus.speed(i),'filled');
     
     % pause 2/10 second: 
     pause(0.05)
end





end



[Kennfeld_Darstellung_Vorne,x_contourf_Vorne,y_contourf_Vorne]= Erstellung_Kennfeld_Darstellung (Motor_Vorne);
[Kennfeld_Darstellung_Hinten,x_contourf_Hinten,y_contourf_Hinten]= Erstellung_Kennfeld_Darstellung (Motor_Hinten);

Dynamic_Plot_2_Motors(Kennfeld_Darstellung_Vorne,x_contourf_Vorne,y_contourf_Vorne,DATA(2,:),DATA(1,:),Kennfeld_Darstellung_Hinten,x_contourf_Hinten,y_contourf_Hinten,DATA(5,:),DATA(4,:),dc)
