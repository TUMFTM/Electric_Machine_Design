%% Animate solution space of NSGA-II optimization

if ~exist('result', 'var')
    error('Result Struct aus Optimierung laden')
end
% if ~exist('Ergebnis', 'var')
%     error('Ergebnis Struct als Referenz aus VSim laden')
% end
% if ~exist('Param', 'var')
%     error('Param Struct als Referenz aus VSim laden')
% end

% % Daten zusammenfassen
% if Param.transmission.ratios(end) < 1
%     Overdrive = 1;
% else
%     Overdrive = 0;
% end
% if Param.transmission.shift_time == 0
%     DSG = 1;
% else
%     DSG = 0;
% end

% Nr. des ausgewählten Individuums
%no_ind = 1;

%ref = [Param.transmission.ratios(1)/Param.transmission.ratios(end), length(Param.transmission.ratios), Overdrive, DSG, Param.final_drive.ratio, Param.engine.M_max, NaN(1,23), Param.TCO.Gesamtkosten/(Param.TCO.Jahresfahrleistung(1)*(Param.vehicle.payload/100000)), (Ergebnis.OUT_summary.signals(7).values(end)*Param.engine.fuel.co2_per_litre*1000 + (-Ergebnis.delta_WPT -Ergebnis.delta_E) * Param.em.fuel.co2_per_kwh) /(Ergebnis.OUT_summary.signals(2).values(length(Ergebnis.OUT_summary.signals(2).values))/1000)/(Param.vehicle.payload/1000), Elastizitaet_auslesen(Param.VSim)];
% EURO VI Referenz Diesel für Hybrid LKW Optimierung
%ref = [Param.transmission.ratios(1)/Param.transmission.ratios(end), length(Param.transmission.ratios), Overdrive, DSG, Param.final_drive.ratio, Param.engine.M_max, NaN(1,23), Param.TCO.Gesamtkosten/(Param.TCO.Jahresfahrleistung(1)*(Param.vehicle.payload/100000)), (Ergebnis.OUT_summary.signals(7).values(end)*Param.engine.fuel.co2_per_litre*1000 + (-Ergebnis.delta_WPT -Ergebnis.delta_E) * Param.em.fuel.co2_per_kwh) /(Ergebnis.OUT_summary.signals(2).values(length(Ergebnis.OUT_summary.signals(2).values))/1000)/(Param.vehicle.payload/1000), Elastizitaet_auslesen(Param.VSim)];
% EURO VI Referenz Diesel für Elektro LKW Optimierung
%ref = [NaN(1,5), Param.final_drive.ratio, Param.transmission.ratios(1)/Param.transmission.ratios(end), length(Param.transmission.ratios), Overdrive, DSG, Param.engine.shift_parameter.n_lo, Param.engine.shift_parameter.n_pref,   NaN(1,2), Param.TCO.Gesamtkosten/(Param.TCO.Jahresfahrleistung(1)*(Param.vehicle.payload/100000)), (Ergebnis.OUT_summary.signals(7).values(end)*Param.engine.fuel.co2_per_litre*1000 + (-Ergebnis.delta_WPT -Ergebnis.delta_E) * Param.em.fuel.co2_per_kwh) /(Ergebnis.OUT_summary.signals(2).values(length(Ergebnis.OUT_summary.signals(2).values))/1000)/(Param.vehicle.payload/1000), Elastizitaet_auslesen(Param.VSim)];
obj_name = {'Kosten in €', 'Zykluseffizienz in %','Maschinenvolumen in m^3'};


% Preallocate
ind = zeros(result.opt.popsize, result.opt.numVar+result.opt.numObj, result.opt.maxGen);
% Results auslesen und einzelne Individuen zusammenfassen
for j=1:result.opt.popsize
    for i = 1:result.opt.maxGen
    ind(j,1:result.opt.numVar+result.opt.numObj,i) = [result.pops(i,j).var, result.pops(i,j).obj];
    end
end

% x_bound = [ 0 10];
% y_bound = [ 0 100];
% z_bound = [ 0 50];
% l_bound = [0 0 0];
% u_bound = [30 100 60];



%% Plotten
fig = figure('Name', 'NSGA-II Animation');
set(gcf,'color','w');
% Größe des Fensters auf Papiergröße anpassen
set(fig, 'Units', 'Pixels');
% DIN A4
set(fig, 'Position', [100, 100, 1000, 1000]);
    
for i=1:result.opt.maxGen
    row_ind = 0;
    %     %-------Optimierung-------------
    %     plot3(ind(:,result.opt.numVar+1,i), ind(:,result.opt.numVar+2,i), ind(:,result.opt.numVar+3,i), 'o','MarkerSize',6, 'color', 'b');
    %     %--------Refernenz------------
    %     plot3(ref(1,result.opt.numVar+1), ref(1,result.opt.numVar+2), ref(1,result.opt.numVar+3),'d', 'MarkerSize', 8, 'LineWidth', 1.5, 'color', 'r');
    clf
    for j = 1:result.opt.numObj
        for k = 1:result.opt.numObj
            % Dont plot diagonal
             if j>k
                subplot(result.opt.numObj ,result.opt.numObj,row_ind+k+1)
                plot(ind(:,result.opt.numVar+j,i), ind(:,result.opt.numVar+k,i), 'o');
                % axis limits
                %axis([l_bound(j) u_bound(j) l_bound(k) u_bound(k)])
                box on
                xlabel(obj_name(j), 'FontSize', 12)
                ylabel(obj_name(k), 'FontSize', 12)
             end

        end
        row_ind = row_ind + k;
    end
    %     % axis limits
    %     axis([x_bound y_bound z_bound]);
    %     view([57, 28]);
    %     grid on
    %     legend({'Population', 'Reference'}, 'FontSize', 9, 'location','best', 'Orientation', 'Vertical');
    drawnow
    M1(i) = getframe(fig);
end


%%
% close all
% [h, w, p] = size(M1(1).cdata);  % use 1st frame to get dimensions
% hf = figure; 
% % resize figure based on frame's w x h, and place at (150, 150)
% set(hf, 'position', [150 150 w h]);
% axis off
% movie(hf,M1);
% %%
% mplay(M1)

%% Lösungsraum plotten

% Welche Generation soll ausgewertet werden
gen = result.opt.maxGen;

% Preallocate
ind = zeros(result.opt.popsize, result.opt.numVar+result.opt.numObj);

% Results auslesen und einzelne Individuen zusammenfassen
for j=1:result.opt.popsize
    ind(j,1:result.opt.numVar+result.opt.numObj) = [result.pops(gen,j).var, result.pops(gen,j).obj];
end

az = linspace(0, 75, 20*10);

fig2 = figure('Name', 'Lösungsraum');
set(gcf,'color','w');
% Größe des Fensters auf Papiergröße anpassen
set(fig2, 'Units', 'Pixels');
% DIN A4
set(fig2, 'Position', [100, 100, 1000, 1000]);

for i = 1:length(az)
clf
%-------Optimierung-------------
hold on;
plot3(ind(:,result.opt.numVar+1), ind(:,result.opt.numVar+2), ind(:,result.opt.numVar+3), 'o','MarkerSize',6, 'color', 'b');
%--------Referenz------------
%plot3(ref(1,result.opt.numVar+1), ref(1,result.opt.numVar+2), ref(1,result.opt.numVar+3),'d', 'MarkerSize', 8, 'LineWidth', 1.5);
%---Ausgewähltes Individuum---
%plot3(ind(no_ind,result.opt.numVar+1), ind(no_ind,result.opt.numVar+2), ind(no_ind,result.opt.numVar+3),'x', 'MarkerSize', 8, 'LineWidth', 1.5);

x1 = xlabel(obj_name(1));
y1 = ylabel(obj_name(2));
z1 = zlabel(obj_name(3)');
set(x1, 'Units', 'Normalized', 'FontSize', 12);
set(y1, 'Units', 'Normalized', 'FontSize', 12, 'Rotation', 0);
set(z1, 'Units', 'Normalized', 'Position', [0, 0.5, 1], 'FontSize', 12, 'Rotation', 90);

hAxis = gca;
hAxis.ZRuler.FirstCrossoverValue  = 8; % Z crossover with X axis
hAxis.ZRuler.SecondCrossoverValue = 45; % Z crossover with Y axis
set(hAxis,  'projection', 'perspective');
% xlim([8 9]);
% ylim([30 45]);

%axislabel_translation_slider;
%legend({'Population', 'Reference', 'Choice'}, 'FontSize', 9, 'location','best', 'Orientation', 'Vertical');
view([az(i), 15]);
grid on;
drawnow
M2(i) = getframe(fig2);
end

%%
% close all
% [h, w, p] = size(M2(1).cdata);  % use 1st frame to get dimensions
% hf = figure; 
% % resize figure based on frame's w x h, and place at (150, 150)
% set(hf, 'position', [150 150 w h]);
% axis off
% movie(hf,M2);
% %%
 mplay(M2)

%% Save movie
writerObj = VideoWriter('NSGA_Animation', 'MPEG-4');
open(writerObj);
% Write out all the frames.

    writeVideo(writerObj, M1);

close(writerObj);
%% Save movie
writerObj = VideoWriter('NSGA_Animation_Loesungsraum', 'MPEG-4');
open(writerObj);
% Write out all the frames.

    writeVideo(writerObj, M2);

close(writerObj);
