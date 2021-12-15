function [MotMap]=CreationMapPlot(map)


B=map>1 ;
MotMap = 1./map.*B + ~B.*map;



end