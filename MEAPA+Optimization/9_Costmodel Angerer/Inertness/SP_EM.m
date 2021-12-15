function [ x_EM, y_EM, z_EM ] = SP_EM( l_EM,l_GTR,b_GTR,akt_az)
%Focus estimation (amounts only, allocation VA/HA in par_MDT).
if akt_az == 1
    x_EM = -0.4*l_GTR;
    y_EM = -(0.5*b_GTR + 0.5*l_EM);
    z_EM = 0;
else %rn must be distinguished left/right: x,y,z = [xyz_left xyz_right], y_axis points to the right wheel, x to the back
    x_EM = [-0.4*l_GTR -0.4*l_GTR];
    y_EM = [-0.5*l_EM  +0.5*l_EM]; 
    z_EM = [0 0];
end

end

