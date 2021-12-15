function [P_EL_MIN] = BatteryEfficiency(P_EL_MIN)
 if P_EL_MIN>=0
     P_EL_MIN=P_EL_MIN/0.97
 else
     P_EL_MIN=P_EL_MIN*0.97
 end
end

