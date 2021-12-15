% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Modell for the design and analysis of PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (svenja.kalt@tum.de), 
%           Jonathan Erhard 
% -------------------------------------------------------------------------


% Electric sheet materials load
list = dir('5_Material/1_ElectricalSheet');
if isempty(list)
    Elektroblech_List = [{'-'}];
else
    Elektroblech_List = {};
    for i = 1:length(list)
        if(length(list(i).name)>4 && strcmp(list(i).name(end-3:end),'.mat'))
            Elektroblech_List = [Elektroblech_List, {list(i).name(1:end-4)}];
        end
    end
    if isempty(Elektroblech_List)
        Elektroblech_List = [{'-'}];
    end
end

% Conductor materials load
list = dir('5_Material/2_Conductor');
if isempty(list)
    Leiter_List = [{'-'}];
else
    Leiter_List = {};
    for i = 1:length(list)
        if(length(list(i).name)>4 && strcmp(list(i).name(end-3:end),'.mat'))
            Leiter_List = [Leiter_List, {list(i).name(1:end-4)}];
        end
    end
    if isempty(Leiter_List)
        Leiter_List = [{'-'}];
    end
end

% Magnet Material load
list = dir('5_Material/3_Magnet');
if isempty(list)
    Magnet_List = [{'-'}];
else
    Magnet_List = {};
    for i = 1:length(list)
        if(length(list(i).name)>4 && strcmp(list(i).name(end-3:end),'.mat'))
            Magnet_List = [Magnet_List, {list(i).name(1:end-4)}];
        end
    end
    if isempty(Magnet_List)
        Magnet_List = [{'-'}];
    end
end

% Machine topology
opt.Maschinenausfuehrung = [{'IPMSM (tangential)'}; {'IPMSM (V-Form)'}; {'IPMSM (eingelassen)'}; {'SPMSM'}];

% Circuit
opt.Schaltung = [{'Stern'}; {'Dreieck'}];

% Coil form Stator
opt.Spulenform_Stator = [{'Runddraht'}];% {'Formspule- oder Stab'}];

% Coil form Rotor
opt.Spulenform_Rotor = [{'-'}];

% Slot form Stator
opt.Nutform_Stator = [{'Trapezform (eckig)'}];% {'Trapezform (rund)'}];

% Slot form Rotor
opt.Nutform_Rotor = [{'-'}];

% Cooling Stator
opt.Kuehlungsart = [{'Wasser (direkt)'}; {'Luft (indirekt)'}];

% Iron material Stator
opt.Stator_Eisenmaterial.String = Elektroblech_List;

% Conductor material Stator
opt.Stator_Leitermaterial.String = Leiter_List;

% Iron material Rotor
opt.Rotor_Eisenmaterial.String = Elektroblech_List;

% Conductor material Rotor
opt.Rotor_Leitermaterial.String = Leiter_List;

% Magnet material Rotor
opt.Rotor_Magnetmaterial.String = Magnet_List;

% Modus winding design
opt.Mode_Wicklung = [{'Klassisch'}; {'Optimierung'}; {'Manuell'}];

% Optimization goal winding
opt.Wicklungstyp = [{'A'}; {'B'}; {'C'}];
