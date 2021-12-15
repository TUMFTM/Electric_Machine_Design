function    [lambda]= calc_Lambda(p);
%% Determination of the relative armature length lambda HERZOG: Script for the lecture: "Entwurf elektrischer Maschinen". SS15. S. 208
    if p==1
        lambda= 1.4;                                         % between 1 and 4
    elseif p>1
        lambda=nthroot(p,3);
    end