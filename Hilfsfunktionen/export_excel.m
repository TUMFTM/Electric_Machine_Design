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

function export_table = export_excel(input)
    var = fieldnames(input);
    
    for i = 1:length(var)
        export_struct{i,1} = input.(var{i});
    end
    
    export_table = table(export_struct,'VariableNames',{'var'},'RowNames',var);
end