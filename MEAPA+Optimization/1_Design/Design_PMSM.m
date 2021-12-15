% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Modell for the design and analysis of PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (svenja.kalt@tum.de), 
%           Jonathan Erhard 
% -------------------------------------------------------------------------


%% Table of contents
% A) Preprocessing
% B) Determination of main dimensions
% B.1) Stator inner diameter
% B.2) Geometric air gap and outer rotor diameter
% B.3) Axial length of active parts
% C) Stator design
% C.1) Winding scheme
% C.2) Wiring
% C.3) Slot shape and magnetic circuit
% D) Rotor design
% D.1) Magnet dimensioning
% D.2) Magnet positioning and magnetic circuit
% E) Recalculation
% E.1) Inductances
% E.2) resistances
% F) Postprocessing
% G) Auxiliary functions

function [Entwurf] = Design_PMSM(handles)
%% A) Preprocessing
% #########################################################################
% #   A) PREPROCESSING                                                    #
% #########################################################################

warning('off','backtrace')

rated = handles.rated;
richt = handles.richt;
opt = handles.opt;


const.mu_0 = pi * 4e-7;

%% B) Main dimensions
% #########################################################################
% #   B) MAIN DIMENSIONS                                                  #
% #########################################################################


rated.M_N = rated.P_N * 60 / (2 * pi * rated.n_N);


if(strcmp(opt.Schaltung,'Dreieck'))
    emag.U_1Str = rated.U_N;
elseif(strcmp(opt.Schaltung,'Stern'))
    emag.U_1Str = rated.U_N / sqrt(3);
else
    error('Ungueltige Eingabe bei Variable "opt.Schaltung"');
end


emag.E_h = emag.U_1Str;

eta_N_data = [ 0.1  0.2  0.4  0.6  0.8  1.0  2.0  4.0  6.0  8.0 10.0 20.0 40.0 60.0 80.0 100.0 200.0 400.0 600.0 800.0 1000.0; ...
              89.7 91.6 93.4 94.4 94.9 95.3 96.3 97.0 97.3 97.5 97.6 98.0 98.3 98.4 98.5  98.6  98.7  98.8  98.9  98.9   98.9; ...
              86.6 88.9 90.8 91.9 92.6 93.1 94.5 95.4 95.9 96.1 96.3 96.8 97.2 97.5 97.7  97.8  98.1  98.3  98.4  98.4  98.4];
emag.eta_N = interp1(eta_N_data(1,:), eta_N_data(2,:),(rated.P_N/1e3),'linear')/100;
if(any(isnan(emag.eta_N)))
    error('Datenbereich der hinterlegten Kurve zur Abschaetzung des Nennwirkungsgrads zu klein. \n(keine Werte fuer den Nennwirkungsgrad eta_N bei einer Nenneistung von P_N = %.2f W hinterlegt)',rated.P_N(1,end));
end
clear eta_N_data


emag.P_el_N = rated.P_N / emag.eta_N;


emag.P_s = rated.P_N / (emag.eta_N * rated.cos_phi_N);


emag.I_1Str = emag.P_s / (rated.m * emag.U_1Str);


if(strcmp(opt.Schaltung,'Dreieck'))
    emag.I_N = emag.I_1Str * sqrt(3);
elseif(strcmp(opt.Schaltung,'Stern'))
    emag.I_N = emag.I_1Str;
else
    error('Ungueltige Eingabe bei Variable "opt.Schaltung"');
end


emag.P_si = emag.P_s * (emag.E_h / emag.U_1Str);


%% B.1) Stator inner diameter
% #########################################################################
% #   B.1) STATOR INNER DIAMETER                                          #
% #########################################################################


C_s_data = [0.5 1.00 2.00 5.00 10.00 20.00 50.00 100.00 200.00 500.00 1000.00; ...
            2.0 2.75 3.11 3.40  3.73  4.00  4.54   4.99   5.49   6.22    6.84];
        
geo.misc.C_s = interp1(C_s_data(1,:), C_s_data(2,:),((emag.P_s/1e3)/(2*rated.p)),'linear');
if(any(isnan(geo.misc.C_s)))
    error('Datenbereich der hinterlegten Kurve zur Abschaetzung des Ausnutzungsfaktors zu klein. \n(keine Werte fuer den Ausnutzungsfaktors C_s bei einer Scheinleistung von P_s = %.2f W hinterlegt)',rated.P_N(1,end));
end
clear C_s_data


geo.misc.C_s = geo.misc.C_s * (0.8 + 0.25 * rated.cos_phi_N);


geo.misc.C = geo.misc.C_s * (emag.E_h / emag.U_1Str);


geo.D_1i = nthroot(((emag.P_s/1e3) * 2*rated.p) / (rated.n_N * geo.misc.C_s * richt.lambda * pi),3);
D_1i_alt = nthroot(((emag.P_si/1e3) * 2*rated.p) / (rated.n_N * geo.misc.C * richt.lambda * pi),3);

if((geo.D_1i-D_1i_alt)>1e-3)
    error('Anomalie in der Berechnung des Bohrungsdurchmessers "D_1i"')
end
clear D_1i_alt


k_s = 3;
v = geo.D_1i * pi * k_s * (rated.n_N/60);
if(v>150)
    warning('maximale Umfangsgeschwindigkeit v_max zu gross, lambda anpassen');
end
clear k_s v

%% B.2) Geometric air gap and outer rotor diameter
% #########################################################################
% #   B.2) GEOMETRIC AIR GAP AND OUTER ROTOR DIAMETER                     #
% #########################################################################


if(rated.p==1)
    geo.delta = (0.2 + 0.01*rated.P_N^0.4);
else
    geo.delta = (0.18 + 0.006*rated.P_N^0.4);
end


if(strcmp(opt.Maschinenausfuehrung,'SPMSM'))
    geo.delta = geo.delta *1.6;
end


if(geo.delta<0.2)
    geo.delta = 0.2;
end


geo.D_2a = geo.D_1i - 2*(geo.delta/1e3);

%% B.3) Axial length of active parts
% #########################################################################
% #   B.3) AXIAL LENGTH OF ACTIVE PARTS                                   #
% #########################################################################


geo.tau_1p = (geo.D_1i * pi) / (2 * rated.p);


geo.l_i = (emag.P_s/1e3) / (geo.misc.C_s * geo.D_1i^2 * rated.n_N);
l_i_alt = geo.tau_1p * richt.lambda;


if((geo.l_i-l_i_alt)>0.001)
    error('Anomalie in der Berechnung der ideellen Laenge "l_i"')
end
clear l_i_alt


if(strcmp(opt.Kuehlungsart,'Luft (indirekt)'))
    if(geo.l_i>0.2)
        geo.misc.n_v = floor(geo.l_i / 0.08);
    else
        geo.misc.n_v = 0;
    end
elseif(strcmp(opt.Kuehlungsart,'Wasser (direkt)'))
    geo.misc.n_v = 0;
else
    error('Ungueltige Eingabe bei Variable "opt.Kuehlungsart"')
end


geo.misc.gamma_v = (1 / (1 + 5 * (geo.delta/1e3 / richt.l_v)));


geo.l = geo.l_i + (geo.misc.n_v * geo.misc.gamma_v * richt.l_v) - 2*(geo.delta/1e3);


geo.l_Fe = geo.l - (geo.misc.n_v * richt.l_v);


geo.V_Bohrung = geo.D_1i^2 * geo.l * (pi/4);

%% C) Stator
% #########################################################################
% #   C) STATOR                                                           #
% #########################################################################


if(isfield(richt, 'B_p'))
    
    emag.A = (geo.misc.C * sqrt(2)) / (pi^2 * richt.B_p * richt.xi_1p) * 60;

    
    if(strcmp(opt.Kuehlungsart,'Luft (indirekt)'))
        if(emag.A<30 || emag.A>120)
            warning('Variable "emag.A" ausserhalb der Grenzen')
        end
    elseif(strcmp(opt.Kuehlungsart,'Wasser (direkt)'))
        if(emag.A<160 || emag.A>300) 
            warning('Variable "emag.A" ausserhalb der Grenzen')
        end
    else
        error('Ungueltige Eingabe bei Variable "opt.Kuehlungsart"')
    end
elseif(isfield(richt, 'A'))

    emag.B_p = (geo.misc.C * sqrt(2)) / (pi^2 * richt.A * richt.xi_1p) * 60;
    
    
    if(emag.B_p<0.75 || emag.B_p>1.05) 
        warning('Variable "emag.B_p" ausserhalb der Grenzen')
    end
    
    richt.B_p = emag.B_p;
    emag.A = richt.A;
else
    error('Kein Wert fuer die Hauptwellen-Amplitude der Luftspaltinduktion "B_p" oder den Strombelag "A" gefunden.')
end


emag.Phi_h = 2/pi * richt.B_p * geo.tau_1p * geo.l_i;


emag.Phi_delta = emag.Phi_h;


emag.Phi_1r_max = emag.Phi_delta / 2;


if(emag.A*richt.S_1<100 || emag.A*richt.S_1>350) 
    warning('Produkt (Strombelag * Stromdichte) ausserhalb der Grenzen')
end


wick.l_1m = 2.0 * (geo.l + 1.3 * geo.tau_1p + (0.03 + 0.02*(emag.U_1Str/1000)));


geo.Nut.h_1r_min = (emag.Phi_1r_max) / (geo.l_Fe * richt.phi_1Fe * richt.B_1r_max);


geo.Nut.h_1r_max = (emag.Phi_1r_max / (2/pi)) / (geo.l_Fe * richt.phi_1Fe * richt.B_1r_max);


geo.Nut.h_1n_min = emag.A / (richt.S_1 * richt.phi_1n * (1-(((richt.B_p)*geo.l_i)/(richt.B_1z_max*richt.phi_1Fe*geo.l_Fe))))/1e3;


geo.Nut.h_1n_max = 1.25 * geo.Nut.h_1n_min;


geo.D_1a_max = geo.D_1i + 2*geo.Nut.h_1n_max + 2*geo.Nut.h_1r_max;

%% C.1) Winding scheme
% #########################################################################
% #   C.1) WINDING SCHEME                                                 #
% #########################################################################


wick.N_1max = floor((geo.D_1i * pi) / richt.tau_1n_min);


N_1 = (1:wick.N_1max)';


q_1 = N_1 ./ (2.*rated.p.*rated.m);


[q_1z, q_1n] = numden(sym(q_1));
q_1z = double(q_1z);
q_1n = double(q_1n);


for i = 1:length(q_1n)
    ggT_q_1n_m(i,1) = ggT_fun(q_1n(i,1),rated.m);
end


tau_1n = (geo.D_1i.*pi) ./ N_1;


y_1Durchmesser = N_1 ./ (2.*rated.p);


wick.table_1all = table(N_1, q_1, q_1z, q_1n, ggT_q_1n_m, tau_1n,  y_1Durchmesser);
clear N_1 q_1 q_1z q_1n ggT_q_1n_m tau_1n y_1Durchmesser


wick.q_min_GL = 2;


wick.N_min_GL = (2 * rated.p * rated.m * wick.q_min_GL);


if(strcmp(opt.Mode_Wicklung,'Klassisch'))
   
    if(wick.N_min_GL>wick.N_1max && ~strcmp(opt.Wicklungstyp,'C'))
        warning('Wicklungstyp A und B nicht moeglich! Wicklungstyp C wird verwendet!');
        opt.Wicklungstyp = 'C';
    end
    
    switch opt.Wicklungstyp
        case 'A'
            Wicklungstypen = {'1SGL'};
        case 'B'
            Wicklungstypen = {'2SGL'};
        case 'C'
            Wicklungstypen = {'2SBL'};
    end
elseif(strcmp(opt.Mode_Wicklung,'Optimierung'))
    Wicklungstypen = {'1SGL', '1SBL', '2SGL', '2SBL'};
elseif(strcmp(opt.Mode_Wicklung,'Manuell'))
    Wicklungstypen = {'1SGL', '2SGL', '2SBL'}; 
    
else
    error('Ungueltige Eingabe bei Variable "opt.Mode_Wicklung"');
end

for i_Wicklung = 1:length(Wicklungstypen)
    
    switch Wicklungstypen{i_Wicklung}
        case {'1SGL', '2SGL'}
            
            W = wick.table_1all(wick.table_1all.q_1n==1,:);
        case {'1SBL', '2SBL'}
            
            W = wick.table_1all(wick.table_1all.q_1n~=1,:);
    end

    
    W.Wicklungstyp_1(:) = Wicklungstypen(i_Wicklung);
    W = [W(:,end) W(:,1:end-1)];

    
    switch Wicklungstypen{i_Wicklung}
        case {'1SGL', '1SBL'}
            W.n_lay = 1 * ones(length(W.N_1),1);
        case {'2SGL', '2SBL'}
            W.n_lay = 2 * ones(length(W.N_1),1);
    end

    
    switch Wicklungstypen{i_Wicklung}
        case '1SBL'
            
            Symmetrie = (~mod((rated.p),W.q_1n) & W.ggT_q_1n_m==1);

            
            W = W(Symmetrie,:);

        case '2SBL'
            
            Symmetrie = (~mod((2*rated.p),W.q_1n) & W.ggT_q_1n_m==1);

           
            W = W(Symmetrie,:);
    end
    clear Symmetrie

    
    switch Wicklungstypen{i_Wicklung}
        case '1SBL'
            W = W(W.q_1<1,:);
    end
    if(isempty(W))
        continue;
    end

    
    switch Wicklungstypen{i_Wicklung}
        case {'1SGL', '2SGL'}
            W = W(W.q_1>=wick.q_min_GL,:);
    end
    if(isempty(W))
        continue;
    end

   
    for i = 1:length(W.N_1)
        count_Sehnung = 0;
        if(W.q_1(i)<1) 
            y_pos = 1;
        else 
            switch Wicklungstypen{i_Wicklung}
                case {'1SGL'}
                    
                    y_pos = W.y_1Durchmesser(i);
                case {'1SBL'}
                    
                case {'2SGL'}
                    y_pos = 1:1:W.y_1Durchmesser(i)-1; 
                case {'2SBL'}
                    y_pos = 1:1:(W.y_1Durchmesser(i)*2);
            end
        end
        y_v_pos = W.y_1Durchmesser(i) - y_pos;
        for j = 1:length(y_v_pos)
            
            if((y_pos(j)/W.y_1Durchmesser(i))>=(2/3) && (y_pos(j)/W.y_1Durchmesser(i))<=(4/3))
                if(count_Sehnung>0)
                    W(end+1,:) = W(i,:);
                    W.y_1(end) = y_pos(j);
                    W.y_1v(end) = y_v_pos(j);
                else
                    W.y_1(i) = y_pos(j);
                    W.y_1v(i) = y_v_pos(j);
                end
                count_Sehnung = count_Sehnung + 1;
            end

        end
        clear count_Sehnung y_pos y_v_pos
    end
    W.Sehnung_1 = W.y_1 ./ W.y_1Durchmesser;

   
    W(W.Sehnung_1==0,:) = [];
    if(isempty(W))
        continue;
    end

   
    W = sortrows(W,{'N_1','Sehnung_1'},{'ascend','ascend'});

    
    b_ns_data = [0.05 0.2 0.3 0.45; 1.5 3 3 5];
    b_ns_fun = polyfit(b_ns_data(1,:),b_ns_data(2,:),1);
    NG.b_1ns = polyval(b_ns_fun,geo.D_1i);
    clear b_ns_data b_ns_fun
    
    if(~strcmp(opt.Mode_Wicklung,'Manuell'))
        for i = 1:length(W.N_1)
            
            g = [0,1:100000,-1:-1:-100000];
            nu_1s = rated.p * (1 + ((2 .* rated.m .* g) ./ W.q_1n(i)));
            nu_1s = sort(nu_1s,2,'ComparisonMethod','abs');
            nu_1 = nu_1s./rated.p;

            
            if(mod(W.q_1n(i),2)) 
                p_u = W.q_1n(i);
                Q_u = 2 * rated.m * W.q_1z(i);
            else 
                p_u = W.q_1n(i) / 2;
                Q_u = rated.m * W.q_1z(i);
            end

            if(mod(Q_u,2)) 
                q_1 = (Q_u + rated.m) / (2 * rated.m);
                q_2 = q_1 - 1;
            else 
                q_1 = Q_u / (2 * rated.m);
                q_2 = q_1;
            end

            alpha_Q = (2 * pi * p_u) / Q_u;

            g_min = 0;
            Y = ((g_min * Q_u) + 1) / p_u;
            while(mod(Y,1))
                g_min = g_min + 1;
                Y = ((g_min * Q_u) + 1) / p_u;
            end

            
            xi_1gr = (sin(nu_1 .* alpha_Q .* Y .* (q_1./2)) - cos(nu_1 .* p_u .* pi .* Y) .* sin(nu_1 .* alpha_Q .* Y .* (q_2./2))) ./ ((q_1 + q_2) .* sin(nu_1 .* alpha_Q .* (Y./2)));

         
            xi_1sp = sin(nu_1 .* (pi./2) .* W.Sehnung_1(i));

           

           
            xi_1n = (sin((nu_1s.*(NG.b_1ns./1e3))./geo.D_1i))./((nu_1s.*(NG.b_1ns./1e3))./geo.D_1i);

            
            xi_1spgr = xi_1sp .* xi_1gr;

            
            xi_1 = xi_1sp .* xi_1gr .* xi_1n;

            
            xi_1p = xi_1(nu_1==1);

           
            var1 = nu_1s(nu_1~=1);
            var2 = xi_1spgr(nu_1~=1);
            W.sigma_1o(i) = sum(((rated.p./var1) .* (var2./xi_1spgr(nu_1==1))).^2);

            
            W.xi_1sp(i) = {xi_1sp(1:15)'};
            W.xi_1gr(i) = {xi_1gr(1:15)'};
            W.xi_1n(i) = {xi_1n(1:15)'};
            W.xi_1spgr(i) = {xi_1spgr(1:15)'};
            W.xi_1(i) = {xi_1(1:15)'};
            W.xi_1p(i) = abs(xi_1p);
            W.nu_1(i) = {nu_1(1:15)'};

            clear g nu_1s nu_1 p_u Q_u q_1 q_2 alpha_Q g_min Y xi_1sp xi_1gr xi_1n xi_1spgr xi_1 xi_1p var1 var2
        end
    end

%% C.2) Wiring
% #########################################################################
% #   C.2) WIRING                                                         #
% #########################################################################

    
    switch Wicklungstypen{i_Wicklung}
        case '1SGL'
            
            a_max = rated.p;

            
            a_pos = 1;
            for a = 2:a_max
                if(~mod((rated.p/a),1))
                    a_pos = [a_pos, a];
                end
            end

        case '2SGL'
           
            a_max = 2 * rated.p;

            
            a_pos = 1;
            for a = 2:a_max
                if(~mod((rated.p/a),1))
                    a_pos = [a_pos, a];
                end
            end
        case {'1SBL', '2SBL'}
            W_copy = W;
            for i = 1:length(W.N_1)
                if(mod(W.q_1n(i),2)) 
                    t = rated.p / W.q_1n(i);

                    
                    a_max = 2*t;

                else 
                    t = 2 * (rated.p / W.q_1n(i));

                    
                    a_max = t;

                end

               
                a_pos = 1;
                for a = 2:a_max
                    if(~mod((t/a),1))
                        a_pos = [a_pos, a];
                    end
                end

                
                count_Parallel = 0;
                for j = 1:length(a_pos)
                    if(count_Parallel>0)
                        W_copy(end+1,:) = W_copy(i,:);
                        W_copy.a_1(end) = a_pos(j);
                    else
                        W_copy.a_1(i) = a_pos(j);
                    end
                    count_Parallel = count_Parallel + 1;
                end

            end
            W = W_copy; 
    end

    switch Wicklungstypen{i_Wicklung}
        case {'1SGL','2SGL'}
            
            for i = 1:length(W.N_1)
                count_Parallel = 0;
                for j = 1:length(a_pos)
                    if(count_Parallel>0)
                        W(end+1,:) = W(i,:);
                        W.a_1(end) = a_pos(j);
                    else
                        W.a_1(i) = a_pos(j);
                    end
                    count_Parallel = count_Parallel + 1;
                end
            end
    end
    clear a_max a_pos a count_Parallel i j t W_copy

    
    W = sortrows(W,{'N_1','Sehnung_1','a_1'},{'ascend','ascend','ascend'});
    
    if(strcmp(opt.Mode_Wicklung,'Manuell') && length(Wicklungstypen)~=i_Wicklung)
    elseif(strcmp(opt.Mode_Wicklung,'Manuell') && length(Wicklungstypen)==i_Wicklung)
        if(i_Wicklung==1)
            wick.table_1pos = W;
        else
            wick.table_1pos = [wick.table_1pos; W];
        end
        clear W saveNG
        
        popup_String = '';
        for i = 1:length(wick.table_1pos.N_1)
            popup_String = [popup_String; {[wick.table_1pos.Wicklungstyp_1{i} ', N_1=' num2str(wick.table_1pos.N_1(i)) ', q_1=' num2str(wick.table_1pos.q_1(i)) ', Sehnung_1=' num2str(wick.table_1pos.Sehnung_1(i)) ', a_1=' num2str(wick.table_1pos.a_1(i))]}];
        end

        handles.popupmenu_Auswahl_Wicklung.String = ['-'; popup_String];
        handles.popupmenu_Auswahl_Wicklung.Value = 1;
        set(handles.uipanel_Auswahl_Wicklung,'Position',[0,0,742,826])
        set(handles.uipanel_Auswahl_Wicklung,'Visible','on')
        set(handles.uipanel_Ergebnisse,'Visible','off')
        handles.text_Wicklungslayout.Visible = 'on';
        handles.edit_Wicklungslayout.Visible = 'on';
        handles.text_y_1v.Visible = 'on';
        handles.edit_y_1v.Visible = 'on';
        handles.pushbutton_Wicklungslayout.Visible = 'on';
        handles.axes_Auswahl_Wicklung.Visible = 'off';
        handles.toolbar.Children(1).Enable = 'off';
        handles.toolbar.Children(2).Enable = 'off';
        handles.toolbar.Children(3).Enable = 'off';
        handles.toolbar.Children(4).Enable = 'off';
        save(['3_Results/1_Misc/Entwurf_',opt.file_id,'_temp.mat'],'wick','rated');

        uiwait(gcf);

        load(['3_Results/1_Misc/Entwurf_',opt.file_id,'_temp.mat']);
        delete(['3_Results/1_Misc/Entwurf_',opt.file_id,'_temp.mat']);

        handles.toolbar.Children(1).Enable = 'on';
        handles.toolbar.Children(2).Enable = 'on';
        handles.toolbar.Children(3).Enable = 'on';
        handles.toolbar.Children(4).Enable = 'on';

        if(handles.popupmenu_Auswahl_Wicklung.Value~=1)
            wick.table_1sel = wick.table_1pos(handles.popupmenu_Auswahl_Wicklung.Value-1,:);
        else
            
            wick.table_1sel.a_1 = 1;
        end
        W = wick.table_1sel;

        set(handles.uipanel_Auswahl_Wicklung,'Position',[1484,0,742,826])
        set(handles.uipanel_Auswahl_Wicklung,'Visible','off')
        set(handles.uipanel_Ergebnisse,'Visible','on')
        clear popup_String
        pause(0.1)

        cla(handles.axes_Auswahl_Wicklung)
        handles.popupmenu_Auswahl_Wicklung.Value = 1;
        handles.popupmenu_Auswahl_Wicklung.String = {'-'};
        handles.edit_y_1v.String = {'-'};

        Wicklungstypen = W.Wicklungstyp_1;
        i_Wicklung = 1 ;

        
        g = [0,1:10000,-1:-1:-10000];
        nu_1s = rated.p * (1 + ((2 .* rated.m .* g) ./ W.q_1n));
        nu_1s = sort(nu_1s,2,'ComparisonMethod','abs');
        nu_1 = nu_1s./rated.p;

       
        h = 1; 
        x=zeros(length(nu_1),1,length(wick.Matrix(h,:)));
        y=zeros(length(nu_1),1,length(wick.Matrix(h,:)));
        z=zeros(length(nu_1),1);
        for i = 1:length(nu_1)
            
            x(i,h,:) = sin(((2.*pi)./W.N_1) .* rated.p .* nu_1(i) .* (abs(wick.Matrix(h,:))-1)) .* sign(wick.Matrix(h,:));
            
            y(i,h,:) = cos(((2.*pi)./W.N_1) .* rated.p .* nu_1(i) .* (abs(wick.Matrix(h,:))-1)) .* sign(wick.Matrix(h,:));
            
            z(i,h) = sqrt((sum(x(i,h,:))).^2 + (sum(y(i,h,:))).^2);
        end

       
        xi_1spgr = (z ./ ((W.N_1*W.n_lay)/rated.m))';

        clear x y z

        
        b_ns_data = [0.05 0.2 0.3 0.45; 1.5 3 3 5];
        b_ns_fun = polyfit(b_ns_data(1,:),b_ns_data(2,:),1);
        NG.b_1ns = polyval(b_ns_fun,geo.D_1i);
        clear b_ns_data b_ns_fun

       
        xi_1n = (sin((nu_1s.*(NG.b_1ns./1e3))./geo.D_1i))./((nu_1s.*(NG.b_1ns./1e3))./geo.D_1i);

        
        xi_1 = xi_1spgr .* xi_1n;

        
        xi_1p = xi_1(nu_1==1);

       
        var1 = nu_1s(nu_1~=1);
        var2 = xi_1spgr(nu_1~=1);
        W.sigma_1o = sum(((rated.p./var1) .* (var2./xi_1spgr(nu_1==1))).^2);

        
        W.xi_1n = {xi_1n(1:15)'};
        W.xi_1spgr = {xi_1spgr(1:15)'};
        W.xi_1 = {xi_1(1:15)'};
        W.xi_1p = abs(xi_1p);
        W.nu_1 = {nu_1(1:15)'};

        clear i h g nu_1s nu_1 xi_1n xi_1spgr xi_1 xi_1p var1 var2
    
    else
    end
    
    if(strcmp(opt.Mode_Wicklung,'Manuell') && length(Wicklungstypen)~=i_Wicklung)
    else
        
        W.w_1Str_opt = (sqrt(2) .* emag.E_h) ./ (2 .* pi .* rated.f_N .* W.xi_1p .* emag.Phi_h);

        
        switch Wicklungstypen{i_Wicklung}
            case {'1SGL', '1SBL'}
                W.z_1n = round((W.w_1Str_opt .* 2 .* W.a_1 .* rated.m) ./ W.N_1);
            case {'2SGL', '2SBL'}
                var = floor((W.w_1Str_opt .* 2 .* W.a_1 .* rated.m) ./ W.N_1);
                W.z_1n = ceil((W.w_1Str_opt .* 2 .* W.a_1 .* rated.m) ./ W.N_1);
                W.z_1n(logical(~mod(var,2))) = var(logical(~mod(var,2)));
                clear var
        end

        
        W(W.z_1n<1,:) = [];
        if(isempty(W))
            continue;
        end

       
        W.w_1Str = (W.N_1 .* W.z_1n) ./ (2 .* W.a_1 .* rated.m);

        
        W.err_1 = W.w_1Str_opt - W.w_1Str;

       
        W.Phi_h = (sqrt(2) .* emag.E_h) ./ (2 .* pi .* rated.f_N .* W.xi_1p .* W.w_1Str);

        
        W.Phi_delta = W.Phi_h;

       
        W.B_p = W.Phi_h ./ ((2./pi) .* geo.tau_1p .* geo.l_i);

       
        W.A = (2 .* W.w_1Str .* rated.m .* emag.I_1Str) ./ (pi .* geo.D_1i.*1e3);

       
        W.S_1 = (emag.A .* richt.S_1) ./ W.A;

%% C.3) Slot shape and magnetic circuit
% #########################################################################
% #   C.3) slot shape and magnetic circuit                                #
% #########################################################################

        
        W.I_1zw = emag.I_1Str ./ W.a_1;

        
        W.A_1L = W.I_1zw ./ W.S_1;

       
        W.d_1L = sqrt((4 .* W.A_1L) ./ pi);

        
        W.V_1Le = ((W.A_1L.*1e-6) .* wick.l_1m .* W.w_1Str .* W.a_1 .* rated.m);
        W.m_1Le = W.V_1Le .* opt.Stator_Leitermaterial.rho_Le;

       
        W.A_1n = (W.z_1n .* W.A_1L) ./ richt.phi_1n;

       
        W.Phi_1r_max = (1.05.*W.Phi_delta) ./ 2;

        
        h_ns_data = [0.05 0.2 0.3 0.45; 0.5 1 1 2];
        h_ns_fun = polyfit(h_ns_data(1,:),h_ns_data(2,:),1);
        NG.h_1ns = polyval(h_ns_fun,geo.D_1i);
        clear h_ns_data h_ns_fun

        
        NG.maxIter = 1e6;
        NG.iter = 0;

        
        saveNG = NG;

        
        if(strcmp(opt.Nutform_Stator,'Trapezform (eckig)'))
            for i = 1:length(W.N_1)
               
                NG = saveNG;

               
                if((0.2 * W.tau_1n(i)*1e3)<NG.b_1ns)
                    NG.b_1n_u = NG.b_1ns;                                    
                else
                    NG.b_1n_u = 0.2 * W.tau_1n(i)*1e3;
                end
                NG.b_1n_o = NG.b_1n_u;
                NG.b_1n_m = (NG.b_1n_o + NG.b_1n_u) / 2;

                
                NG.d_1iso = 0.3;

                
                NG.alpha_1nk = pi/6;

                
                NG.h_1k = tan(NG.alpha_1nk) * 0.5 * (NG.b_1n_u - NG.b_1ns);

                
                NG.h_1nk = NG.h_1ns + NG.h_1k + 0.5;

               
                NG.h_1n = NG.h_1nk + 1;

               
                NG.b_1z_o = (geo.D_1i*1e3 + 2*NG.h_1n) * pi / W.N_1(i) - NG.b_1n_o;

                
                NG.b_1z_u = (geo.D_1i*1e3 + 2*NG.h_1ns + 2*NG.h_1k) * pi / W.N_1(i) - NG.b_1n_u;

                
                NG.b_1z_m = (NG.b_1z_o + NG.b_1z_u) / 2;

               
                NG.h_1r = (geo.D_1a_max - geo.D_1i)/2*1e3 - NG.h_1n;

                
                W.A_1n_tat(i) = 0.5 * (NG.b_1n_u + NG.b_1n_o) * (NG.h_1n - NG.h_1nk);

                while(W.A_1n_tat(i)<W.A_1n(i) && NG.iter<NG.maxIter)
                    
                    W.B_1z_o(i) = (W.B_p(i) * W.tau_1n(i) * geo.l_i) / (NG.b_1z_o/1e3 * richt.phi_1Fe * geo.l_Fe);

                   
                    W.B_1z_u(i) = (W.B_p(i) * W.tau_1n(i) * geo.l_i) / (NG.b_1z_u/1e3 * richt.phi_1Fe * geo.l_Fe);

                   
                    W.B_1z_m(i) = (W.B_p(i) * W.tau_1n(i) * geo.l_i) / (NG.b_1z_m/1e3 * richt.phi_1Fe * geo.l_Fe);

                    
                    W.B_1r(i) = W.Phi_1r_max(i) / (NG.h_1r*1e-3 * richt.phi_1Fe * geo.l_Fe);

                    
                    if(W.B_1z_o(i) < richt.B_1z_max || W.B_1z_u(i) < richt.B_1z_max)
                        if(W.B_1z_o(i) < W.B_1z_u(i))
                            NG.b_1n_o = NG.b_1n_o + 0.1;
                        else
                            NG.b_1n_u = NG.b_1n_u + 0.1;
                        end
                    elseif(W.B_1r(i) < W.B_1z_o(i) && W.B_1r(i) < W.B_1z_u(i))
                        NG.h_1n = NG.h_1n + 0.1;
                    elseif(W.B_1z_o(i) < W.B_1z_u(i))
                        NG.b_1n_o = NG.b_1n_o + 0.1;
                    else
                        NG.b_1n_u = NG.b_1n_u + 0.1;
                    end

                    
                    NG.b_1n_m = (NG.b_1n_o + NG.b_1n_u) / 2;
                    NG.b_1z_o = (geo.D_1i*1e3 + 2*NG.h_1n) * pi / W.N_1(i) - NG.b_1n_o;
                    NG.b_1z_u = (geo.D_1i*1e3 + 2*NG.h_1ns + 2*NG.h_1k) * pi / W.N_1(i) - NG.b_1n_u;
                    NG.b_1z_m = (NG.b_1z_o + NG.b_1z_u) / 2;
                    NG.h_1r = (geo.D_1a_max - geo.D_1i)/2*1e3 - NG.h_1n/1e3;
                    NG.h_1k = tan(NG.alpha_1nk) * 0.5 * (NG.b_1n_u - NG.b_1ns);
                    NG.h_1nk = NG.h_1ns + NG.h_1k + 0.5;

                   
                    if NG.b_1z_o<=0 || NG.b_1z_u<=0 || NG.b_1z_m<=0 || NG.h_1r<=0
                        error('Geometrie ist entartet');
                    end

                    
                    W.A_1n_tat(i) = 0.5 * (NG.b_1n_u + NG.b_1n_o) * (NG.h_1n - NG.h_1nk);

                    NG.iter = NG.iter + 1;
                end

                
                if(NG.iter>NG.maxIter)
                    error('Kontrolle Nutraumbilanz notwendig!')
                end

                
                if(W.B_1r(i)<richt.B_1r_max)
                    NG.h_1r = W.Phi_1r_max(i) / (geo.l_Fe * richt.phi_1Fe * richt.B_1r_max) * 1e3;
                end

                
                W.B_1r(i) = W.Phi_1r_max(i) / (NG.h_1r*1e-3 * richt.phi_1Fe * geo.l_Fe);

                
                NG.h_1r_rel = NG.h_1r*1e-3 / geo.tau_1p;

                
                W.D_1a(i) = geo.D_1i + 2 * (NG.h_1r*1e-3 + NG.h_1n*1e-3);

               
                W.Nut_1(i) = NG;

                clear NG
            end
        else
            error('Ungueltige Eingabe bei Variable "opt.Nutform_Stator"');
        end
    end

    if(i_Wicklung==1)
        wick.table_1pos = W;
    else
        wick.table_1pos = [wick.table_1pos; W];
    end
    clear W saveNG

end

clear Wicklungstypen i_Wicklung i j

if(~isfield(wick,'table_1pos'))
    error('Keine gueltige Wicklung gefunden"');
end

if(strcmp(opt.Mode_Wicklung,'Klassisch'))
    
    switch opt.Wicklungstyp
        case 'A'
            wick.table_1sel = wick.table_1pos(wick.table_1pos.Sehnung_1==1,:);
        case {'B', 'C'}
            
            Sehnung_opt = 5/6;
            wick.table_1sel = wick.table_1pos(wick.table_1pos.Sehnung_1==Sehnung_opt,:);
            if(isempty(wick.table_1sel))
                Sehnung_opt = 2/3;
                wick.table_1sel = wick.table_1pos(wick.table_1pos.Sehnung_1==Sehnung_opt,:);
                if(isempty(wick.table_1sel))
                    error('Keine gueltige Wicklung gefunden"');
                end
            end
    end
    
    
    if(max(wick.table_1pos.q_1)<1) 
        wick.table_1sel = wick.table_1sel(:,wick.table_1sel.N_1<=(3*rated.p) & wick.table_1sel.N_1(1,:)>=(1.5*rated.p));

        
        [~,idx] = min(abs((1/rated.m) - wick.table_1sel.q_1));
        wick.table_1sel = wick.table_1sel(idx,:);
    else
        
        wick.table_1sel = wick.table_1sel(max(wick.table_1sel.N_1)==wick.table_1sel.N_1,:);
        
        
        [~, idx] = min(wick.table_1sel.err_1);
        wick.table_1sel = wick.table_1sel(idx,:);
    end
    var1 = table2struct(wick.table_1sel);
    var2 = [fieldnames(wick)' fieldnames(var1)'; struct2cell(wick)' struct2cell(var1)'];
    wick = struct(var2{:});

    clear Sehnung_opt idx var1 var2
    
    wick = TingleyAlg(rated, wick);
    
elseif(strcmp(opt.Mode_Wicklung,'Optimierung'))
    
    del=0;
    for i = 1:length(wick.table_1pos.N_1)
        for j = 1:length(wick.table_1pos.N_1)
            if(i==j)
            else
            if(wick.table_1pos.sigma_1o(i)>=wick.table_1pos.sigma_1o(j) && wick.table_1pos.xi_1p(i)<=wick.table_1pos.xi_1p(j) && abs(wick.table_1pos.err_1(i))>=abs(wick.table_1pos.err_1(j)) && wick.table_1pos.D_1a(i)>=wick.table_1pos.D_1a(j) && wick.table_1pos.m_1Le(i)>=wick.table_1pos.m_1Le(j))
                if(wick.table_1pos.sigma_1o(i)==wick.table_1pos.sigma_1o(j) && wick.table_1pos.xi_1p(i)==wick.table_1pos.xi_1p(j) && abs(wick.table_1pos.err_1(i))==abs(wick.table_1pos.err_1(j)) && wick.table_1pos.D_1a(i)==wick.table_1pos.D_1a(j) && wick.table_1pos.m_1Le(i)==wick.table_1pos.m_1Le(j))
                else
                    del(i) = i;
                    break;
                end
            end
            end
        end
    end
    del = unique(del);
    del = del(del~=0);
    del = flipud(del');
    wick.table_1sel = wick.table_1pos;
    for i = 1:length(del)
        wick.table_1sel(del(i),:) = [];
    end    
    
   
    wick.table_1sel = wick.table_1sel(abs(wick.table_1sel.err_1)<=((wick.table_1sel.w_1Str_opt*1.1)-wick.table_1sel.w_1Str_opt),:);
    
    if(length(wick.table_1sel.N_1)>7)
       
        if(max(wick.table_1sel.xi_1p)>=0.9)
            wick.table_1sel = wick.table_1sel(wick.table_1sel.xi_1p>=0.9,:);
        end
    end
    
    if(length(wick.table_1sel.N_1)>7)
       
        if(min(wick.table_1sel.sigma_1o)<=0.1 && max(wick.table_1sel.q_1)>1)
            wick.table_1sel = wick.table_1sel(wick.table_1sel.sigma_1o<=0.1,:);
        end
    end
    
   
    
    set(handles.figure1, 'currentaxes', handles.axes_Auswahl_Wicklung);
    cla(handles.axes_Auswahl_Wicklung)
    for i = 1:length(wick.table_1sel.N_1)
        P(i,:) = [wick.table_1sel.xi_1p(i) wick.table_1sel.sigma_1o(i)*1e2 wick.table_1sel.err_1(i) wick.table_1sel.D_1a(i)*1e3 wick.table_1sel.m_1Le(i)];
    end
    P_labels = [{'xi_{1p}'};{'sigma_{1o} in %'};{'err_1 (abs)'};{'D_{1a} in mm'};{'m_{1Le} in kg'}];
    axes_interval = 5;
	axes_precision = 3;
    spider_plot(P, P_labels, axes_interval, axes_precision,...
                'Marker', 'o',...
                'LineWidth', 2,...
                'MarkerSize', 5);
            
    % Title properties
    title(handles.axes_Auswahl_Wicklung,'Wicklungen','interpreter','latex','FontSize', 20)
    
    % Legend properties
    lgd_String = '';
    for i = 1:length(wick.table_1sel.N_1)
        lgd_String = [lgd_String; {[wick.table_1sel.Wicklungstyp_1{i} ', N_1=' num2str(wick.table_1sel.N_1(i)) ', q_1=' num2str(wick.table_1sel.q_1(i)) ', Sehnung_1=' num2str(wick.table_1sel.Sehnung_1(i)) ', a_1=' num2str(wick.table_1sel.a_1(i))]}];
    end
    lgd = legend(lgd_String, 'Location', 'southoutside');
    if(length(wick.table_1sel.N_1)>10)
        lgd.NumColumns = 2;
    end
    
  
    handles.popupmenu_Auswahl_Wicklung.String = ['-'; lgd_String];
    handles.popupmenu_Auswahl_Wicklung.Value = 1;
    set(handles.uipanel_Auswahl_Wicklung,'Position',[0,0,742,826])
    set(handles.uipanel_Auswahl_Wicklung,'Visible','on')
    set(handles.uipanel_Ergebnisse,'Visible','off')
    handles.text_Wicklungslayout.Visible = 'off';
    handles.edit_Wicklungslayout.Visible = 'off';
    handles.pushbutton_Wicklungslayout.Visible = 'off';
    handles.text_y_1v.Visible = 'off';
    handles.edit_y_1v.Visible = 'off';
    handles.axes_Auswahl_Wicklung.Visible = 'off';
    handles.toolbar.Children(1).Enable = 'off';
    handles.toolbar.Children(2).Enable = 'off';
    handles.toolbar.Children(3).Enable = 'off';
    handles.toolbar.Children(4).Enable = 'off';

    uiwait(gcf);
    
    wick.table_1sel = wick.table_1sel(handles.popupmenu_Auswahl_Wicklung.Value-1,:);
    var1 = table2struct(wick.table_1sel);
    var2 = [fieldnames(wick)' fieldnames(var1)'; struct2cell(wick)' struct2cell(var1)'];
    wick = struct(var2{:});
    
    set(handles.uipanel_Auswahl_Wicklung,'Position',[1484,0,742,826])
    set(handles.uipanel_Auswahl_Wicklung,'Visible','off')
    set(handles.uipanel_Ergebnisse,'Visible','on')
    handles.toolbar.Children(1).Enable = 'on';
    handles.toolbar.Children(2).Enable = 'on';
    handles.toolbar.Children(3).Enable = 'on';
    handles.toolbar.Children(4).Enable = 'on';
    
    cla(handles.axes_Auswahl_Wicklung,'reset')
    handles.popupmenu_Auswahl_Wicklung.Value = 1;
    handles.popupmenu_Auswahl_Wicklung.String = {'-'};
    
    clear i j del P P_labels axes_interval axes_precision lgd_String var1 var2
    
    wick = TingleyAlg(rated, wick);
    
elseif(strcmp(opt.Mode_Wicklung,'Manuell'))
    var1 = table2struct(wick.table_1pos);
    temp = wick.Matrix_lay;
    wick = rmfield(wick,'Matrix_lay');
    var2 = [fieldnames(wick)' fieldnames(var1)'; struct2cell(wick)' struct2cell(var1)'];
    wick = struct(var2{:});
    wick.Matrix_lay = temp;
    
    clear var1 var2 temp
else
    error('Ungueltige Eingabe bei Variable "opt.Mode_Wicklung"');
end

% Data handling
emag.w_1Str_opt = wick.w_1Str_opt;
wick = rmfield(wick,'w_1Str_opt');
emag.w_1Str = wick.w_1Str;
wick = rmfield(wick,'w_1Str');
emag.Phi_h = wick.Phi_h;
wick = rmfield(wick,'Phi_h');
emag.Phi_delta = wick.Phi_delta;
wick = rmfield(wick,'Phi_delta');
emag.B_p = wick.B_p;
wick = rmfield(wick,'B_p');
emag.A = wick.A;
wick = rmfield(wick,'A');
emag.I_1zw = wick.I_1zw;
wick = rmfield(wick,'I_1zw');
emag.Phi_1r_max = wick.Phi_1r_max;
wick = rmfield(wick,'Phi_1r_max');
emag.B_1z_o = wick.B_1z_o;
wick = rmfield(wick,'B_1z_o');
emag.B_1z_u = wick.B_1z_u;
wick = rmfield(wick,'B_1z_u');
emag.B_1z_m = wick.B_1z_m;
wick = rmfield(wick,'B_1z_m');
emag.B_1r = wick.B_1r;
wick = rmfield(wick,'B_1r');
geo.Nut.A_1n = wick.A_1n;
wick = rmfield(wick,'A_1n');
geo.Nut.A_1n_tat = wick.A_1n_tat;
wick = rmfield(wick,'A_1n_tat');
geo.D_1a = wick.D_1a;
wick = rmfield(wick,'D_1a');
geo.Nut_1 = wick.Nut_1;
wick = rmfield(wick,'Nut_1');

var2 = [fieldnames(geo.Nut_1)' fieldnames(geo.Nut)'; struct2cell(geo.Nut_1)' struct2cell(geo.Nut)'];
geo.Nut_1 = struct(var2{:});
geo = rmfield(geo,'Nut');
clear var2

%% D) Rotor
% #########################################################################
% #   D) ROTOR                                                            #
% #########################################################################


geo.D_2i = geo.D_2a / 3;


emag.Phi_2r_max = emag.Phi_delta / 2;


geo.tau_2p = (geo.D_2a * pi) / (2 * rated.p);


emag.misc.gamma_1 = 1 / (1 + 5 * geo.delta / (geo.Nut_1.b_1ns));
emag.misc.k_1c = wick.tau_1n / (wick.tau_1n - emag.misc.gamma_1* geo.Nut_1.b_1ns*1e-3);


emag.misc.k_c = emag.misc.k_1c;


geo.delta_i = emag.misc.k_c * geo.delta;

%% D.1) Magnet dimensioning
% #########################################################################
% #   D.1) MAGNET DIMENSIONING                                            #
% #########################################################################


geo.l_PM = geo.l_Fe;


if(strcmp(opt.Maschinenausfuehrung,'IPMSM (V-Form)'))
geo.anzahl_PM = 4*rated.p;
else
geo.anzahl_PM = 2*rated.p; 
end

if(strcmp(opt.Maschinenausfuehrung,'SPMSM'))
    
    geo.h_2r = (geo.D_2a-geo.D_2i)/2 *1e3;
        
    
    emag.B_2r = emag.Phi_2r_max / (geo.h_2r*1e-3 * richt.phi_1Fe * geo.l_Fe);
    
   
    emag.misc.k_so = 1.1; 
    
    
    geo.alpha_i = 2/3; 
    
    
    geo.b_PM = geo.alpha_i * geo.tau_2p *1e3;
    
   
    emag.misc.k_sigma = 1.0; 
    
    
    emag.misc.k_b = 1; 
    
    
    emag.misc.k_h = 2; 
    
   
    emag.misc.k_vg = pi/2; 
    
   
    emag.misc.k_h_q = 1; 
    
    
    emag.misc.k_sd = 1.1; 
    emag.misc.k_sq = 1.1; 
    
    
    geo.h_PM = ((emag.B_p * 2 * opt.Rotor_Magnetmaterial.mu_r * emag.misc.k_so * (geo.delta_i*1e-3)) / ...
        ((opt.Rotor_Magnetmaterial.B_r - (emag.B_p * (geo.alpha_i / emag.misc.k_sigma) * ((pi * geo.D_1i) / (2 * rated.p * emag.misc.k_b * geo.b_PM*1e-3)))) * emag.misc.k_h)) * 1e3;
    
   
    geo.D_2a = geo.D_2a - 2*geo.h_PM*1e-3;
    
    
    geo.tau_2p = (geo.D_2a * pi) / (2 * rated.p);
    
    
    geo.delta_i = geo.delta_i + (geo.h_PM/opt.Rotor_Magnetmaterial.mu_r);
else
   
    emag.Phi_PM = emag.Phi_delta;
    
    
    geo.b_PM = (emag.Phi_PM / (emag.B_p * geo.l_PM)) * 1e3;
    
    if(strcmp(opt.Maschinenausfuehrung,'IPMSM (eingelassen)'))
       
        geo.h_2r = ((geo.D_2a)/2) * 1e3;
    
    elseif(strcmp(opt.Maschinenausfuehrung,'IPMSM (tangential)'))
        
        geo.Abstand_PM_Rotoroberflaeche = 3;
        
        
        geo.h_1_PM = ((geo.D_2a)/2 *1e3) - 0.5*sqrt((geo.D_2a*1e3)^2 - geo.b_PM^2);
        
        
        geo.h_2r = ((geo.D_2a)/2 *1e3) - (geo.h_1_PM + geo.Abstand_PM_Rotoroberflaeche);
        
    elseif(strcmp(opt.Maschinenausfuehrung,'IPMSM (V-Form)'))
        
        geo.alpha_PM = 120 * ((2*pi)/360);
        
        
        geo.Abstand_PM_unten = 3;
        
       
        geo.Abstand_PM_Rotoroberflaeche = 3;
        
       
        geo.h_1_PM = ((geo.D_2a)/2 *1e3) - 0.5*sqrt((geo.D_2a*1e3)^2 - ((geo.b_PM * sin(geo.alpha_PM/2)) + geo.Abstand_PM_unten)^2);
        geo.h_2_PM = (geo.b_PM/2) * cos(geo.alpha_PM/2);
        
        
        geo.h_2r = ((geo.D_2a)/2 *1e3) - (geo.h_1_PM + geo.Abstand_PM_Rotoroberflaeche + geo.h_2_PM);
    else
        error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"');
    end
    
    
    emag.B_2r = emag.Phi_2r_max / (geo.h_2r*1e-3 * richt.phi_1Fe * geo.l_Fe);
    
    
    emag.V_delta = (1/const.mu_0) * emag.B_p * emag.misc.k_c * geo.delta*1e-3;

    
    emag.B_1z_o = (emag.B_p * wick.tau_1n * geo.l_i) / (geo.Nut_1.b_1z_o*1e-3 * richt.phi_1Fe * geo.l_Fe);
    emag.B_1z_m = (emag.B_p * wick.tau_1n * geo.l_i) / (geo.Nut_1.b_1z_m*1e-3 * richt.phi_1Fe * geo.l_Fe);
    emag.B_1z_u = (emag.B_p * wick.tau_1n * geo.l_i) / (geo.Nut_1.b_1z_u*1e-3 * richt.phi_1Fe * geo.l_Fe);
    
    
    emag.H_1z_o = interp1(opt.Stator_Eisenmaterial.B, opt.Stator_Eisenmaterial.H, emag.B_1z_o);
    emag.H_1z_m = interp1(opt.Stator_Eisenmaterial.B, opt.Stator_Eisenmaterial.H, emag.B_1z_m);
    emag.H_1z_u = interp1(opt.Stator_Eisenmaterial.B, opt.Stator_Eisenmaterial.H, emag.B_1z_u);
    if(any(isnan(emag.H_1z_o)))
        error('Datenbereich der hinterlegten BH-Kurve (Elektroband) zu klein. Bitte anderes Eisenmaterial waehlen. \n(keine Werte fuer die Zahnfeldstaerke (Stator) H_1z_o bei einer Zahninduktion (Stator) von B_1z_o = %.2f T hinterlegt)',emag.B_1z_o(1,end));
    end
    if(any(isnan(emag.H_1z_m)))
        error('Datenbereich der hinterlegten BH-Kurve (Elektroband) zu klein. Bitte anderes Eisenmaterial waehlen. \n(keine Werte fuer die Zahnfeldstaerke (Stator) H_1z_m bei einer Zahninduktion (Stator) von B_1z_m = %.2f T hinterlegt)',emag.B_1z_m(1,end));
    end
    if(any(isnan(emag.H_1z_u)))
        error('Datenbereich der hinterlegten BH-Kurve (Elektroband) zu klein. Bitte anderes Eisenmaterial waehlen. \n(keine Werte fuer die Zahnfeldstaerke (Stator) H_1z_u bei einer Zahninduktion (Stator) von B_1z_u = %.2f T hinterlegt)',emag.B_1z_u(1,end));
    end
    
    
    emag.V_1z = (1/6) * (emag.H_1z_o + 4*emag.H_1z_m + emag.H_1z_u) * geo.Nut_1.h_1n*1e-3;

    
    emag.V_z = emag.V_1z;

    
    emag.H_1r = interp1(opt.Stator_Eisenmaterial.B, opt.Stator_Eisenmaterial.H, emag.B_1r);
    if(any(isnan(emag.H_1r)))
        error('Datenbereich der hinterlegten BH-Kurve (Elektroband) zu klein. Bitte anderes Eisenmaterial waehlen. \n(keine Werte fuer die Rueckenfeldstaerke (Stator) H_1r bei einer Rueckeninduktion (Stator) von B_1r = %.2f T hinterlegt)',emag.B_1r(1,end));
    end

    
    emag.H_2r = interp1(opt.Rotor_Eisenmaterial.B, opt.Rotor_Eisenmaterial.H, emag.B_2r);
    if(any(isnan(emag.H_2r)))
        error('Datenbereich der hinterlegten BH-Kurve (Elektroband) zu klein. Bitte anderes Eisenmaterial waehlen. \n(keine Werte fuer die Rueckenfeldstaerke (Rotor) H_2r bei einer Rueckeninduktion (Rotor) von B_2r = %.2f T hinterlegt)',emag.B_2r(1,end));
    end

    
    geo.tau_1r = (geo.D_1i + 2*geo.Nut_1.h_1n*1e-3)* pi / (2*rated.p);

    
    geo.tau_2r = (geo.D_2a - (2*geo.h_2r*1e-3 + geo.D_2i))* pi / (2*rated.p);

    
    C_r_data = [0:0.1:2.0; ...
                0.72 0.72 0.72 0.72 0.72 0.72 0.71 0.70 0.68 0.63 0.57 0.49 0.40 0.31 0.25 0.20 0.17 0.15 0.14 0.13 0.12];
    
    emag.misc.C_1r = interp1(C_r_data(1,:), C_r_data(2,:),emag.B_1r,'linear');
    emag.misc.C_2r = interp1(C_r_data(1,:), C_r_data(2,:),emag.B_2r,'linear');
    if(any(isnan(emag.misc.C_1r)))
        error('Datenbereich der hinterlegten Kurve zur Abschaetzung des Rueckenreduktionsfaktors zu klein. \n(keine Werte fuer den Rueckenreduktionsfaktor C_1r bei einem Quotient von h_1r/tau_1r, = %.2f hinterlegt)',geo.Nut_1.h_1r*1e-3/geo.tau_1r);
    end
    if(any(isnan(emag.misc.C_2r)))
        error('Datenbereich der hinterlegten Kurve zur Abschaetzung des Rueckenreduktionsfaktors zu klein. \n(keine Werte fuer den Rueckenreduktionsfaktor C_2r bei einem Quotient von h_2r/tau_2r, = %.2f hinterlegt)',geo.h_2r*1e-3/geo.tau_2r);
    end

   
    emag.V_1r = emag.misc.C_1r * emag.H_1r * (geo.tau_1r/2);

   
    emag.V_2r = emag.misc.C_2r * emag.H_2r * (geo.tau_2r/2);

    
    emag.V_r = emag.V_1r + emag.V_2r;
    
   
    geo.h_PM = ((emag.V_delta + emag.V_1z + emag.V_r) / ...
        ((opt.Rotor_Magnetmaterial.B_r - emag.B_p) /  (opt.Rotor_Magnetmaterial.mu_r * const.mu_0))) * 1e3;
    
    if(strcmp(opt.Maschinenausfuehrung,'IPMSM (eingelassen)'))
        
        geo.h_2r = (((geo.D_2a)/2) * 1e3) - geo.h_PM;
        
    elseif(strcmp(opt.Maschinenausfuehrung,'IPMSM (tangential)'))
        
        geo.h_2r = ((geo.D_2a)/2 *1e3) - (geo.h_1_PM + geo.Abstand_PM_Rotoroberflaeche + geo.h_PM);
        
    elseif(strcmp(opt.Maschinenausfuehrung,'IPMSM (V-Form)'))
        
        geo.h_3_PM = 0.5 * geo.h_PM * sqrt(4-cos(geo.alpha_PM/2)^2);
        
        
        geo.h_2r = ((geo.D_2a)/2 *1e3) - (geo.h_1_PM + geo.Abstand_PM_Rotoroberflaeche + geo.h_2_PM + geo.h_3_PM);
    end
    
    
    emag.misc.k_so = 1.1; 
    
   
    alpha_i = 0.95; 
    
    
    emag.misc.k_sigma = 0.9; 
    
   
    emag.misc.k_b = 2;
    
   
    emag.misc.k_h = 1; 
    
    emag.misc.k_vg = 4/pi; 
    
    
    emag.misc.k_h_q = 0;
    
    
    emag.misc.k_sd = 1.1; 
    emag.misc.k_sq = 1.7;
end

%% E) Recalculation of the magnetic circuit
% #########################################################################
% #   E) RECALCULATION OF THE MAGNETIC CIRCUIT                            #
% #########################################################################


if(geo.Nut_1.h_1r_rel<0.13 || geo.Nut_1.h_1r_rel>0.28)
    warning('h_1r_rel ausserhalb der Grenzen')
end


geo.A_1r = (pi/4) * (geo.D_1a^2 - (geo.D_1i + (2*geo.Nut_1.h_1n*1e-3))^2);

geo.A_1z = (pi/4) * (geo.D_1a^2 - geo.D_1i^2) - geo.A_1r - wick.N_1*(geo.Nut_1.A_1n_tat + geo.Nut_1.b_1ns*geo.Nut_1.h_1ns + geo.Nut_1.b_1n_u*(geo.Nut_1.h_1nk-geo.Nut_1.h_1k-geo.Nut_1.h_1ns) + 0.5*(geo.Nut_1.b_1ns*geo.Nut_1.b_1n_u)*geo.Nut_1.h_1k)/1000000;

geo.A_2r = (pi/4) * (geo.D_2a^2 - geo.D_2i^2);


geo.Vo_1r = geo.A_1r * geo.l_Fe * richt.phi_1Fe; 

geo.Vo_1z = geo.A_1z * geo.l_Fe * richt.phi_1Fe;

geo.Vo_2r = geo.A_2r * geo.l_Fe * richt.phi_2Fe;

%% E.1) Inductance
% #########################################################################
% #   E.1) INDUCTANCE                                                     #
% #########################################################################


emag.L_1h = (rated.m/2) * (const.mu_0/(geo.delta_i*1e-3)) * (2/pi) * geo.tau_1p * geo.l_i * (4/pi) * ((emag.w_1Str*wick.xi_1p)^2/(2*rated.p));


b_1ns_delta = geo.Nut_1.b_1ns / geo.delta;
if(b_1ns_delta<3)
    lambda_z_data = [1.15,1,0.8,0.62,0.5,0.39,0.3,0.26,0.19,0.13,0.1,0.07,0.03,0;...
                    0,0.125,0.3125,0.5,0.75,1,1.3125,1.5,1.75,2,2.25,2.5,2.75,3];
    lambda_z_fun = polyfit(lambda_z_data(2,:),lambda_z_data(1,:),4);
    emag.lambda_1z = polyval(lambda_z_fun,b_1ns_delta);
else
    lambda_z_data = [0,-0.03,-0.05,-0.075,-0.1,-0.11,-0.12,-0.13,-0.14,-0.15,-0.16,-0.17,-0.18;...
                    3,4,5,6,7,8,9,10,11,12,13,14,16];
    lambda_z_fun = polyfit(lambda_z_data(2,:),lambda_z_data(1,:),4);
    emag.lambda_1z = polyval(lambda_z_fun,b_1ns_delta);
end
clear b_1ns_delta b_2ns_delta lambda_z_data lambda_z_fun


geo.Nut_1.h_1ue = geo.Nut_1.h_1nk - geo.Nut_1.h_1ns - geo.Nut_1.h_1k;
geo.Nut_1.h_1l = geo.Nut_1.h_1n - geo.Nut_1.h_1nk - geo.Nut_1.d_1iso;
geo.Nut_1.b_1n = geo.Nut_1.b_1n_u; %(b_n_o + b_n_u) / 2;
geo.Nut_1.b_1k = (geo.Nut_1.b_1n + geo.Nut_1.b_1ns) / 2;

switch wick.Wicklungstyp_1
    case {'1SGL', '1SBL'}
       
        emag.lambda_1nz = (geo.Nut_1.h_1l/(3 * geo.Nut_1.b_1n)) + (geo.Nut_1.h_1ue/geo.Nut_1.b_1n) + (geo.Nut_1.h_1k/geo.Nut_1.b_1k) + (geo.Nut_1.h_1ns/geo.Nut_1.b_1ns) + emag.lambda_1z; 

    case {'2SGL', '2SBL'}
        
        k1_data =  [0.82,0.89,0.95,1; 2/3,0.8,0.9,1];
        k1_fun = polyfit(k1_data(2,:),k1_data(1,:),1);
        emag.misc.k_1 = polyval(k1_fun,wick.Sehnung_1);

        k2_data =  [0.75,0.85,0.925,1; 2/3,0.8,0.9,1];
        k2_fun = polyfit(k2_data(2,:),k2_data(1,:),1);
        emag.misc.k_2 = polyval(k2_fun,wick.Sehnung_1);

        
        geo.Nut_1.d = geo.Nut_1.d_1iso;

       
        emag.lambda_1nz = emag.misc.k_1 * (geo.Nut_1.h_1l/(3 * geo.Nut_1.b_1n)) + emag.misc.k_2 * ((geo.Nut_1.h_1ue/geo.Nut_1.b_1n) + (geo.Nut_1.h_1k/geo.Nut_1.b_1k) + (geo.Nut_1.h_1ns/geo.Nut_1.b_1ns) + emag.lambda_1z) + (geo.Nut_1.d/(4 * geo.Nut_1.b_1n)); 
        
        clear k1_data k1_fun k2_data k2_fun
    otherwise
        error('Ungueltige Eingabe bei Variable Wicklung');
end


emag.L_1sigma_nz = 2 * const.mu_0 * geo.l_i * (emag.w_1Str^2/rated.p) * (emag.lambda_1nz/wick.q_1);


wick.l_1w = 0.5 * (wick.l_1m - 2*geo.l);

switch wick.Wicklungstyp_1
    case {'1SGL', '1SBL'}
        
        emag.lambda_1ws = 0.3;
        
        
        emag.lambda_1w = emag.lambda_1ws * (wick.l_1w / (geo.l_i));
        
    case {'2SGL', '2SBL'}
        
        emag.lambda_1ws = 0.25;
        
        
        emag.lambda_1w = emag.lambda_1ws * (wick.l_1w / (2 * geo.l_i));
    otherwise
        error('Ungueltige Eingabe bei Variable Wicklung');
end


emag.L_1sigma_w = 2 * const.mu_0 * geo.l_i * (emag.w_1Str^2/rated.p) * emag.lambda_1w;


emag.L_1sigma_o = wick.sigma_1o * emag.L_1h;


emag.xi_schr_p = 1;
emag.sigma_schr = 1-emag.xi_schr_p^2;


emag.L_sigma_schr = emag.sigma_schr * emag.L_1h;


emag.L_1sigma = emag.L_1sigma_nz + emag.L_1sigma_w + emag.L_1sigma_o + emag.L_sigma_schr;

if(strcmp(opt.Maschinenausfuehrung,'SPMSM'))
    
    emag.misc.k_ad = (geo.delta_i) / (geo.delta_i + (emag.misc.k_h / 2) * (geo.h_PM / opt.Rotor_Magnetmaterial.mu_r));
    emag.misc.k_aq = (geo.delta_i) / (geo.delta_i + (emag.misc.k_h_q) * (geo.h_PM / opt.Rotor_Magnetmaterial.mu_r));

   
    emag.L_hd = emag.L_1h * (emag.misc.k_ad / emag.misc.k_sd) *1.5;
    emag.L_hq = emag.L_1h * (emag.misc.k_aq / emag.misc.k_sq) *1.5;

    
    emag.L_d = emag.L_hd + emag.L_1sigma;
    emag.L_q = emag.L_hq + emag.L_1sigma;

   
    emag.misc.k_phi = (emag.L_1h+emag.L_1sigma)/emag.L_1h;
    
    
    emag.Psi_PM = rated.M_N / (3/sqrt(2) * emag.I_N * rated.p);
else
   
    emag.L_d = 0;
    emag.L_q = 0;
    emag.Psi_PM = 0;
    
    
    emag.misc.k_ad = (geo.delta_i) / (geo.delta_i + (emag.misc.k_h / 2) * (geo.h_PM / opt.Rotor_Magnetmaterial.mu_r));
    emag.misc.k_aq = (geo.delta_i) / (geo.delta_i + (emag.misc.k_h_q) * (geo.h_PM / opt.Rotor_Magnetmaterial.mu_r));

    
    emag.L_hd = emag.L_1h * (emag.misc.k_ad / emag.misc.k_sd) *1.5;
    emag.L_hq = emag.L_1h * (emag.misc.k_aq / emag.misc.k_sq) *1.5;
    
    
    emag.L_d = emag.L_hd + emag.L_1sigma;
    emag.L_q = emag.L_hq + emag.L_1sigma;
    
    
    emag.misc.k_phi = (emag.L_1h+emag.L_1sigma)/emag.L_1h;
    
    
    emag.Psi_PM = emag.misc.k_phi * emag.L_hd * emag.I_N * sqrt(2);
    %Psi_PM = xi_sz_haupt * w_Str * k_vg * alpha_i * richt.B_delta * (D_i/prim.p) * l_e;
end

%% E.2) Resistance
% #########################################################################
% #   E.2) RESISTANCE                                                     #
% #########################################################################


emag.misc.rho_1 = opt.Stator_Leitermaterial.rho_20 * (1 + opt.Stator_Leitermaterial.alpha * (opt.theta_1 - 20));    


emag.misc.kappa_1 = 1/emag.misc.rho_1;


emag.R_1 = (emag.w_1Str * wick.l_1m) / (wick.a_1 * emag.misc.kappa_1 * wick.A_1L);

%% F) Postprocessing
% #########################################################################
% #   F) POSTPROCESSING                                                   #
% #########################################################################

warning('on','backtrace')


geo = x_y_Koordinaten(rated, geo, wick, opt);


opt.Locked = 1;

const = orderfields(const);
emag = orderfields(emag);
emag.misc = orderfields(emag.misc);
geo = orderfields(geo);
geo.Nut_1 = orderfields(geo.Nut_1);
geo.misc = orderfields(geo.misc);
wick = orderfields(wick);

opt = orderfields(opt);
rated = orderfields(rated);
richt = orderfields(richt);


Entwurf.Konstanten = const;
Entwurf.EMAG = emag;
Entwurf.Geometrie = geo;
Entwurf.Wicklung = wick;

Entwurf.Bemessungswerte = rated;
Entwurf.Richtwerte = richt;
Entwurf.Optionen = opt;

end

%% G) Help functions
% #########################################################################
% #   G) HELP FUNCTIONS                                                   #
% #########################################################################

function x = ggT_fun(a,b)
    while b ~= 0 
        z = mod(a,b); 
        a = b; 
        b = z; 
    end 
    x = a;
end


function wick = TingleyAlg(rated, wick)

    
    M1 = zeros(1,wick.q_1n*wick.N_1);
    for i = 1:wick.N_1
        M1(1,(wick.q_1n*(i-1))+1) = i;
    end
    
    n_z = 2*rated.p;

    
    n_s = wick.q_1n * wick.N_1/(2*rated.p);

    
    M2 = reshape(M1,n_s,n_z)';

    
    vz_odd = ones(n_z,1);
    vz_odd(2:2:end,1) = -1;
    vz_even = ones(n_z,1);
    vz_even(1:2:end,1) = -1;
    for i = 1:rated.m
        if(mod(i,2)) %odd
            M3(:,(i-1)*(n_s/rated.m)+1:i*(n_s/rated.m)) = M2(:,(i-1)*(n_s/rated.m)+1:i*(n_s/rated.m)) .* vz_odd;
        else %even
            M3(:,(i-1)*(n_s/rated.m)+1:i*(n_s/rated.m)) = M2(:,(i-1)*(n_s/rated.m)+1:i*(n_s/rated.m)) .* vz_even;
        end
    end

    
    M4 = M3(M3~=0);
    M4 = reshape(M4,[],rated.m);
    M5 = num2cell(M4,1)';
    
    
    switch wick.Wicklungstyp_1
        case {'1SGL', '1SBL'}
            M7 = M4';
            wick.Matrix_lay = M5;
        case {'2SGL', '2SBL'}
            for i = 1:rated.m
                M6(:,i) = sign(M4(:,i)).*(abs(M4(:,i))+wick.y_1).*(-1);
                M7(i,:) = [M4(:,i); M6(:,i)];
            end
            M8 = num2cell(M6,1)';
            M9 = [M5 M8];
            wick.Matrix_lay = M9;
    end
    M7(abs(M7)>wick.N_1) = M7(abs(M7)>wick.N_1)-sign(M7(abs(M7)>wick.N_1))*wick.N_1;
    wick.Matrix = M7;
    
    clear M1 M2 M3 M4 M5 n_z n_s vz_odd vz_even i
    
end


function geo = x_y_Koordinaten(rated, geo, wick, opt)
    
    theta = linspace(0,2*pi,4e2);
    geo.Stator_aussen_x = (geo.D_1a*1e3)/2 * cos(theta) + 0;
    geo.Stator_aussen_y = (geo.D_1a*1e3)/2 * sin(theta) + 0;

    geo.Stator_aussen_x = geo.Stator_aussen_x';
    geo.Stator_aussen_y = geo.Stator_aussen_y';
    
    
    if(strcmp(opt.Nutform_Stator,'Trapezform (eckig)'))
        Segmenthoehe = ((geo.D_1i*1e3)/2) * (1 - sqrt(1 - (geo.Nut_1.b_1ns/(geo.D_1i*1e3))^2));
        x1_Nut = geo.Nut_1.b_1ns / 2;
        x2_Nut = geo.Nut_1.b_1n_u / 2;
        x3_Nut = geo.Nut_1.b_1n_o / 2;
        y1_Nut = (geo.D_1i*1e3)/2 - Segmenthoehe;
        y2_Nut = geo.Nut_1.h_1ns + (geo.D_1i*1e3)/2 - Segmenthoehe;
        y3_Nut = geo.Nut_1.h_1ns + tan(geo.Nut_1.alpha_1nk) * ((geo.Nut_1.b_1n_u-geo.Nut_1.b_1ns) / 2) + (geo.D_1i*1e3)/2 - Segmenthoehe;
        y4_Nut = geo.Nut_1.h_1n + (geo.D_1i*1e3)/2 - Segmenthoehe;
        Nut = [x1_Nut,x1_Nut,x2_Nut,x3_Nut,-x3_Nut,-x2_Nut,-x1_Nut,-x1_Nut; y1_Nut,y2_Nut,y3_Nut,y4_Nut,y4_Nut,y3_Nut,y2_Nut,y1_Nut]';

        
        x1_Fuellung = x2_Nut - geo.Nut_1.d_1iso;
        x2_Fuellung = x3_Nut - geo.Nut_1.d_1iso;
        y1_Fuellung = geo.Nut_1.h_1nk + (geo.D_1i*1e3)/2 - Segmenthoehe;
        y2_Fuellung = geo.Nut_1.h_1n - geo.Nut_1.d_1iso + (geo.D_1i*1e3)/2 - Segmenthoehe;
        Fuellung = [x1_Fuellung,x2_Fuellung,-x2_Fuellung,-x1_Fuellung; y1_Fuellung,y2_Fuellung,y2_Fuellung,y1_Fuellung]';
    else
        error('Ungueltige Eingabe bei Variable "opt.Nutform_Stator"');
    end
    
    
    theta = linspace((asin(geo.Nut_1.b_1ns/(geo.D_1i*1e3))+(pi/2)),(((2*pi)/wick.N_1)-(asin(geo.Nut_1.b_1ns/(geo.D_1i*1e3)))+(pi/2)),1e1);
    Kreissegment = [(geo.D_1i*1e3)/2 * cos(theta); (geo.D_1i*1e3)/2 * sin(theta)]';
    
    
    geo.Stator_innen_x = [];
    geo.Stator_innen_y = [];
    geo.Fuellung_Stator_x = [];
    geo.Fuellung_Stator_y = [];
    for i = 1:wick.N_1
        theta = -((2*(i-1))/wick.N_1)*pi;
        M = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        
        
        Nut_trans = Nut * M;
        
        
        Fuellung_trans = Fuellung * M;
        
        
        Kreissegment_trans = Kreissegment * M;

        
        geo.Stator_innen_x = [geo.Stator_innen_x; Nut_trans(2:end-1,1); Kreissegment_trans(:,1)];
        geo.Stator_innen_y = [geo.Stator_innen_y; Nut_trans(2:end-1,2); Kreissegment_trans(:,2)];
        geo.Fuellung_Stator_x = [geo.Fuellung_Stator_x; Fuellung_trans(:,1)];
        geo.Fuellung_Stator_y = [geo.Fuellung_Stator_y; Fuellung_trans(:,2)];
        geo.Fuellung_Stator_x(end+1,1) = Fuellung_trans(1,1);
        geo.Fuellung_Stator_y(end+1,1) = Fuellung_trans(1,2);
    end
    geo.Stator_innen_x(end+1,1) = geo.Stator_innen_x(1,1);
    geo.Stator_innen_y(end+1,1) = geo.Stator_innen_y(1,1);
    
    % Magnet
    if(strcmp(opt.Maschinenausfuehrung,'SPMSM'))
        % Rotor 
        theta = [0:((2*geo.b_PM*1e-3)/(geo.D_2a + 2*geo.h_PM*1e-3))/1e2:2*pi, 0];%linspace(0,2*pi,1e3);
        geo.Rotor_aussen_x = (geo.D_2a*1e3)/2 * cos(theta) + 0;
        geo.Rotor_aussen_y = (geo.D_2a*1e3)/2 * sin(theta) + 0;

        geo.Rotor_aussen_x = geo.Rotor_aussen_x';
        geo.Rotor_aussen_y = geo.Rotor_aussen_y';
        
        theta = ((2*geo.b_PM*1e-3)/(geo.D_2a + 2*geo.h_PM*1e-3))/1e2:((2*geo.b_PM*1e-3)/(geo.D_2a + 2*geo.h_PM*1e-3))/1e2:2*pi;
        geo.Magnet_Rotor_x = [];
        geo.Magnet_Rotor_y = [];
        geo.points_Magnet = 1;
        for i = 1:(2*rated.p)
            
            theta_1 = theta(1,(theta-(((2*geo.b_PM*1e-3)/(geo.D_2a + 2*geo.h_PM*1e-3))*i + (((2*geo.tau_2p)/geo.D_2a)-((2*geo.b_PM*1e-3)/(geo.D_2a + 2*geo.h_PM*1e-3)))*(i-1)))<=1e-2) + (pi/2);
            theta((theta-((2*geo.tau_2p)/geo.D_2a)*i)<=0) = [];
            Kreissegment_aussen = [((geo.D_2a*1e3)/2 + geo.h_PM) * cos(theta_1); ((geo.D_2a*1e3)/2 + geo.h_PM) * sin(theta_1)]';
            Kreissegment_innen = [((geo.D_2a*1e3)/2) * cos(theta_1); ((geo.D_2a*1e3)/2) * sin(theta_1)]';
            Kreissegment = [Kreissegment_innen(1,:); Kreissegment_aussen; flipud(Kreissegment_innen)];
            
            geo.points_Magnet(i+1) = geo.points_Magnet(i) + length(Kreissegment(:,1));
            
            geo.Magnet_Rotor_x = [geo.Magnet_Rotor_x; Kreissegment(:,1)];
            geo.Magnet_Rotor_y = [geo.Magnet_Rotor_y; Kreissegment(:,2)];
        end
    elseif(strcmp(opt.Maschinenausfuehrung,'IPMSM (eingelassen)'))
        % Rotor 
        theta = [0:((2*geo.b_PM*1e-3)/geo.D_2a)/1e2:2*pi, 0];%linspace(0,2*pi,1e3);
        geo.Rotor_aussen_x = (geo.D_2a*1e3)/2 * cos(theta) + 0;
        geo.Rotor_aussen_y = (geo.D_2a*1e3)/2 * sin(theta) + 0;

        geo.Rotor_aussen_x = geo.Rotor_aussen_x';
        geo.Rotor_aussen_y = geo.Rotor_aussen_y';
        
        theta = ((2*geo.b_PM*1e-3)/geo.D_2a)/1e2:((2*geo.b_PM*1e-3)/geo.D_2a)/1e2:2*pi;
        geo.Magnet_Rotor_x = [];
        geo.Magnet_Rotor_y = [];
        geo.points_Magnet = 1;
        for i = 1:(2*rated.p)
           
            theta_1 = theta(1,(theta-(((2*geo.b_PM*1e-3)/geo.D_2a)*i + (((2*geo.tau_2p)/geo.D_2a)-((2*geo.b_PM*1e-3)/geo.D_2a))*(i-1)))<=1e-2);
            theta((theta-((2*geo.tau_2p)/geo.D_2a)*i)<=0) = [];
            Kreissegment_innen = [((geo.D_2a*1e3)/2 - geo.h_PM) * cos(theta_1); ((geo.D_2a*1e3)/2 - geo.h_PM) * sin(theta_1)]';
            Kreissegment_aussen = [((geo.D_2a*1e3)/2) * cos(theta_1); ((geo.D_2a*1e3)/2) * sin(theta_1)]';
            Kreissegment = [Kreissegment_innen(1,:); Kreissegment_aussen; flipud(Kreissegment_innen)];
            
            geo.points_Magnet(i+1) = geo.points_Magnet(i) + length(Kreissegment(:,1));
            
            geo.Magnet_Rotor_x = [geo.Magnet_Rotor_x; Kreissegment(:,1)];
            geo.Magnet_Rotor_y = [geo.Magnet_Rotor_y; Kreissegment(:,2)];
        end
    elseif(strcmp(opt.Maschinenausfuehrung,'IPMSM (tangential)'))
        % Rotor
        theta = linspace(0,2*pi,4e2);
        geo.Rotor_aussen_x = (geo.D_2a*1e3)/2 * cos(theta) + 0;
        geo.Rotor_aussen_y = (geo.D_2a*1e3)/2 * sin(theta) + 0;

        geo.Rotor_aussen_x = geo.Rotor_aussen_x';
        geo.Rotor_aussen_y = geo.Rotor_aussen_y';
        
        x1_Magnet = geo.b_PM/2;
        x2_Magnet = -geo.b_PM/2;
        y1_Magnet = geo.h_2r;
        y2_Magnet = geo.h_2r + geo.h_PM;
        Magnet = [x1_Magnet,x1_Magnet,x2_Magnet,x2_Magnet; y1_Magnet,y2_Magnet,y2_Magnet,y1_Magnet]';
        
        geo.Magnet_Rotor_x = [];
        geo.Magnet_Rotor_y = [];
        for i = 1:(2*rated.p)
            theta = -((2*(i-1))/(2*rated.p))*pi;
            M = [cos(theta) -sin(theta); sin(theta) cos(theta)];

            % Magnet
            Magnet_trans = Magnet * M;

            
            geo.Magnet_Rotor_x = [geo.Magnet_Rotor_x; Magnet_trans(:,1)];
            geo.Magnet_Rotor_y = [geo.Magnet_Rotor_y; Magnet_trans(:,2)];
            geo.Magnet_Rotor_x(end+1,1) = Magnet_trans(1,1);
            geo.Magnet_Rotor_y(end+1,1) = Magnet_trans(1,2);
        end
        
    elseif(strcmp(opt.Maschinenausfuehrung,'IPMSM (V-Form)'))
        % Rotor
        theta = linspace(0,2*pi,4e2);
        geo.Rotor_aussen_x = (geo.D_2a*1e3)/2 * cos(theta) + 0;
        geo.Rotor_aussen_y = (geo.D_2a*1e3)/2 * sin(theta) + 0;

        geo.Rotor_aussen_x = geo.Rotor_aussen_x';
        geo.Rotor_aussen_y = geo.Rotor_aussen_y';
        
        x1_Magnet = 0;
        x2_Magnet = -geo.b_PM/2;
        y1_Magnet = 0;
        y2_Magnet = -geo.h_PM;
        Magnet = [x1_Magnet,x2_Magnet,x2_Magnet,x1_Magnet,x1_Magnet; y1_Magnet,y1_Magnet,y2_Magnet,y2_Magnet,y1_Magnet]';
        
        theta = pi/2 - geo.alpha_PM/2;
        M = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        
        Magnet = Magnet * M;
        
        MovePM_y = ((geo.D_2a)/2 *1e3) - (geo.h_1_PM + geo.Abstand_PM_Rotoroberflaeche + geo.h_2_PM);
        MovePM_x = -(geo.Abstand_PM_unten/2);
        MovePM = [MovePM_x MovePM_x MovePM_x MovePM_x MovePM_x; MovePM_y MovePM_y MovePM_y MovePM_y MovePM_y]';
        Magnet = Magnet + MovePM;
        Magnet = [Magnet; [Magnet(:,1).*-1  Magnet(:,2)]];
        
        geo.Magnet_Rotor_x = [];
        geo.Magnet_Rotor_y = [];
        for i = 1:(2*rated.p)
            theta = -((2*(i-1))/(2*rated.p))*pi;
            M = [cos(theta) -sin(theta); sin(theta) cos(theta)];

            % Magnet
            Magnet_trans = Magnet * M;

            
            geo.Magnet_Rotor_x = [geo.Magnet_Rotor_x; Magnet_trans(:,1)];
            geo.Magnet_Rotor_y = [geo.Magnet_Rotor_y; Magnet_trans(:,2)];
        end
    else
        error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"')
    end

    % Rotor 
    theta = linspace(0,2*pi,4e2);
    geo.Rotor_innen_x = (geo.D_2i*1e3)/2 * cos(theta) + 0;
    geo.Rotor_innen_y = (geo.D_2i*1e3)/2 * sin(theta) + 0;

    geo.Rotor_innen_x = geo.Rotor_innen_x';
    geo.Rotor_innen_y = geo.Rotor_innen_y';
end