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

%% Plot
hold(handles.axes_plot,'off');
handles.axes_plot.Visible = 'on';
plot(handles.axes_plot,kenn.omega_k_vec/prim.p/(2*pi)*60,kenn.M_max_vec);
hold(handles.axes_plot,'on');
colorbar(handles.axes_plot,'EastOutside');
grid(handles.axes_plot,'on');
handles.axes_plot.Layer = 'top';
xlabel(handles.axes_plot,'Rot. Speed n in $\frac{U}{min}$','interpreter','latex','FontSize', 15);
ylabel(handles.axes_plot,'Torque M in Nm','interpreter','latex','FontSize', 15);