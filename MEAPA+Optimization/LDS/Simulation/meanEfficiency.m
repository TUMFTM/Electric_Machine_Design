function [meanEffMix, effLoadPointsMix] = meanEfficiency(vehicle, Tn)
%this function calculates the mean efficiency of all load points

    for i = 1:length(vehicle.MOTOR)
        if ~isempty(vehicle.MOTOR{i})
            %intObj = griddedInterpolant({fliplr(vehicle.MOTOR{i}.eff_T_axis),vehicle.MOTOR{i}.eff_n_axis},vehicle.MOTOR{i}.eff); %create interpolation object
            intObj = griddedInterpolant({vehicle.MOTOR{i}.eff_T_axis,vehicle.MOTOR{i}.eff_n_axis},vehicle.MOTOR{i}.eff); %create interpolation object
            
            effMap = (vehicle.MOTOR{i}.eff); %write the eff map from each motor in this new value

            effLoadPoints(i,:) = intObj(Tn(i*2-1,:),Tn(i*2,:)*60); %interpolate eff. for all load points using the interp. object

            effRecu(i,:) = effLoadPoints.*(effLoadPoints>=1); %Recuperation efficiencies (eff. between 0 and 1)
            effRecu(i,:) = 1./effRecu(i,:);
            effRecu(i,effRecu==inf) = 0; 

            effAccc(i,:)  = effLoadPoints.*(effLoadPoints<1); %Acceleration efficiency

            effLoadPointsMix(i,:) = effRecu(i,:)+effAccc(i,:); %this is the mixed efficiency when both recu and acc are between 0 and 1;

            meanEff(i) = mean(effLoadPoints); %mean efficiency of all load points

            meanEffMix(i) = mean(effLoadPointsMix); %mean efficiency of all load points with mixed efficiency
        end
    end
    
end