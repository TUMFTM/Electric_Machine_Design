% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Model for the design of a permanent magnet excited synchronous machine and
% subsequent efficiency map calculation
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (kalt@ftm.mw.tum.de)
%           Jonathan Erhard
%           Prof. Markus Lienkamp
% -------------------------------------------------------------------------

function x = ggT_fun(a,b)
    while b ~= 0 
        z = mod(a,b); 
        a = b; 
        b = z; 
    end 
    x = a;
end