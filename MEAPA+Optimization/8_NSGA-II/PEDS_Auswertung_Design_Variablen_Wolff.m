% Skript zur Auswertung der Optimierung

if ~exist('result', 'var')
    error('Result Struct aus Optimierung laden')
end

% Nr. des ausgewählten Individuums

%no_ind = 1; % X in graphik zu markieren

% Limit für Plots
limit_l = result.opt.lb;
limit_u = result.opt.ub;
d_var_name = {'Nennleistung (kW)', 'Nenndrehzahl (1/min)', 'Spannung (V)', 'Polpaarzahl (-)', 'Maschinentyp (-)', 'Kuehlungsart (-)', 'Magnetanordnung (-)', 'Betteriekapazität (kWh)', 'Getriebeuebersetzung (-)', 'Anzahl Maschinen (-)'}; %X
obj_name = {'Kosten in €', 'Zykluseffizienz in %', 'Maschinenvolumen in m^3'}; %y-Namen

% Welche Generation soll ausgewertet werden
gen = result.opt.maxGen;

% Preallocate
ind = zeros(result.opt.popsize, result.opt.numVar+result.opt.numObj);

% Results auslesen und einzelne Individuen zusammenfassen
for j=1:result.opt.popsize
    ind(j,1:result.opt.numVar+result.opt.numObj) = [result.pops(gen,j).var, result.pops(gen,j).obj];
end

% Designvariablen und Reihenfolge für Optimerung
sort_ind = [1, 2, 3, 4, 5, 6, 7, 8 , 9, 10]; 
%sort_ind = [1:result.opt.numVar]; % wenn andere Reihenfolge gewünscht, einfach 3,2,4,1

% Seitenanzahl
if length(sort_ind) > 6
    no_pages = 2;
    no_plots = [6; length(sort_ind) - 6];
    k = 0;
else
    no_pages = 1;
    no_plots = length(sort_ind);
    k = 0;
end


%% A4 bzw. US Letter Darstellung Plotten

for ii = 1:no_pages
    
    figure('Name', sprintf('Auswertung Optimierung MEAPA %i / %i', ii, no_pages));
    
    % Größe des Fensters auf Papiergröße anpassen
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [20,20,21.59, 27.94 * no_plots(ii) / 6]); % A4 Plot
    movegui('center')
%     set(gcf, 'Position', [1,2,19, 27.94 * no_plots(ii) / 6]); % US Letter
    for j = 1:result.opt.numObj:no_plots(ii) * result.opt.numObj
          
        for i = 0:result.opt.numObj-1
            subplot(no_plots(ii) ,result.opt.numObj,j+i)
            set(gca, 'Units', 'centimeters');
            hold on
            
            if different_color
                for jj = 1:result.opt.popsize
                    if ind(jj,5) == 1 %PMSM
                        plot(ind(jj,result.opt.numVar+1+i), ind(jj,sort_ind(floor(j/7)+1+k)),'.', 'MarkerSize', 20, 'color', 'b');
                    elseif ind(jj,5) == 2 %ASM
                        plot(ind(jj,result.opt.numVar+1+i), ind(jj,sort_ind(floor(j/7)+1+k)),'.', 'MarkerSize', 20, 'color', 'g');
                    else
                        plot(ind(jj,result.opt.numVar+1+i), ind(jj,sort_ind(floor(j/7)+1+k)),'.', 'MarkerSize', 20, 'color', 'r');
                    end
                end 
            else
                plot(ind(:,result.opt.numVar+1+i), ind(:,sort_ind(floor(j/result.opt.numObj)+1+k)), '.', 'MarkerSize', 20); %Hier Größe der Punkte einstellen
            end
            % Ausgewähltes Individuum
            %plot(ind(no_ind,result.opt.numVar+1+i), ind(no_ind,sort_ind(floor(j/3)+1+k)),'x', 'MarkerSize', 8, 'LineWidth', 1.5); %hier wird das rote X geplottet
            
            hold off

            % Beschriftung hinzufügen
            xlabel(obj_name(i+1),'FontName','Times New Roman','FontSize',10)
            y = ylabel(d_var_name{sort_ind(floor(j/result.opt.numObj)+1+k)},'FontName','Times New Roman','FontSize',10);
            set(gca,'FontName','Times New Roman','FontSize',10);
            %ylim([limit_l(sort_ind(floor(j/3)+1+k)) limit_u(sort_ind(floor(j/3)+1+k))]);
            
            % Plot um 0,75cm nach oben verschieben, um Platz für Legende zu
            % haben
            y_temp = get(gca, 'Position');
            %set(gca, 'Position', [y_temp(1), y_temp(2)+0.75, 4, 2.5]);
            
            % Überlappen von ylabel mit Label und Ziffern verhindern
            set(y, 'Units', 'Normalized', 'Position', [-0.175, 0.5, 0]);
        end

    end
    
    %% Darstellung anpassen
    
    % Legende manuell mittig, unten platzieren
%     hL = legend({'Population of HYB1 in the 150th Generation', 'EuroVI Diesel (P_I_C_E=324kW)', 'Optimium Individual of HYB1'}, 'FontSize', 9, 'Orientation', 'Horizontal');
%     newPosition = [0.25 0.05 0.5 0];
%     newUnits = 'centimeters';
%     set(hL,'Position', newPosition,'Units', newUnits);
            
    k = no_plots(ii);
end

%% 3D-Plot
figure;
plot3(ind(:,result.opt.numVar+1), ind(:,result.opt.numVar+2),ind(:,result.opt.numVar+3),'.', 'MarkerSize', 20);
grid on
xlabel('Kosten in €')
ylabel('Zykluseffizienz in %')
zlabel('Maschinenvolumen in m^3')
set(gca,'FontName','Times New Roman','FontSize',10);

%% Unnötige Variablen löschen
clearvars -except result Ergebnis Param