function Dynamic_Plot_1_Motor(Kennfeld_Darstellung,x_contourf, y_contourf,Drehzahl,Moment, Zyklus)

 [TumColors] = tumColors(  );

subplot(2,1,1)
contourf(x_contourf, y_contourf,Kennfeld_Darstellung)
title('Motor 1')

hold on

subplot(2,1,2);
plot (Zyklus.speed)
title('Zyklus')

hold on

     subplot(2,1,1)
     s1=scatter(Drehzahl(1),Moment(1), 'filled');

     subplot(2,1,2);
     s3=scatter(1, Zyklus.speed(1),'filled');
     

for i = 2:length(Drehzahl) 
    
    s1.CData=TumColors.secondaryDarkGrey;
    s3.CData=TumColors.primaryBlue;
    
     s1=scatter(subplot(2,1,1),Drehzahl(i),Moment(i), 'filled');
  
     
     s3=scatter(subplot(2,1,2),i, Zyklus.speed(i),'filled');
     
     % pause 2/10 second: 
     pause(0.05)
end