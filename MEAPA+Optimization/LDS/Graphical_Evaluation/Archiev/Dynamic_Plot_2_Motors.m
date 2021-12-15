function Dynamic_Plot_2_Motors(Kennfeld_Darstellung_1,x_contourf_1, y_contourf_1,Drehzahl_1,Moment_1,Kennfeld_Darstellung_2,x_contourf_2, y_contourf_2,Drehzahl_2,Moment_2, Zyklus)

 [TumColors] = tumColors(  );

subplot(2,2,1)
contourf(x_contourf_1, y_contourf_1,Kennfeld_Darstellung_1)
title('Motor 1')

hold on

subplot(2,2,2)
contourf(x_contourf_2, y_contourf_2,Kennfeld_Darstellung_2)
title('Motor 2')

hold on

subplot(2,2,3:4);
plot (Zyklus.speed)
title('Zyklus')

hold on

     subplot(2,2,1)
     s1=scatter(Drehzahl_1(1),Moment_1(1), 'filled');
     
     subplot(2,2,2)
     s2=scatter(Drehzahl_2(1),Moment_2(1), 'filled');
     
     subplot(2,2,3:4);
     s3=scatter(1, Zyklus.speed(1),'filled');
     

for i = 2:length(Drehzahl_1) 
    
    s1.CData=TumColors.secondaryDarkGrey;
    s2.CData=TumColors.secondaryDarkGrey;
    s3.CData=TumColors.primaryBlue;
    
     s1=scatter(subplot(2,2,1),Drehzahl_1(i),Moment_1(i), 'filled');
     
     s2=scatter(subplot(2,2,2),Drehzahl_2(i),Moment_2(i), 'filled');
     
     s3=scatter(subplot(2,2,3:4),i, Zyklus.speed(i),'filled');
     
     % pause 2/10 second: 
     pause(0.05)
end