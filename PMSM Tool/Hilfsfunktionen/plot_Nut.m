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

%% groove form
% Trapezoidal slot coordinates
% x-zero point lies exactly in the center of the slot
% y-zero point lies in groove opening

x1 = b_ns / 2;
x2 = b_n_u / 2;
x3 = b_n_o / 2;
y1 = 0;
y2 = h_ns;
y3 = h_ns + tan(alpha_nk) * (b_n_u-b_ns) / 2;
y4 = h_n;

% Plot groove form
opt.axes_Animate_Nut.Visible = 'on';
plot(opt.axes_Animate_Nut,[x1,x1,x2,x3,-x3,-x2,-x1,-x1],[y1,y2,y3,y4,y4,y3,y2,y1],'k')

%% Slot fill factor
% Fill Coordinates
x1 = x2 - d_iso;
x2 = x3 - d_iso;
y1 = h_nk;
y2 = h_n - d_iso;

% Plot slot form
hold(opt.axes_Animate_Nut, 'on')
fill(opt.axes_Animate_Nut,[x1,x2,-x2,-x1],[y1,y2,y2,y1],'r')

%% Adjustment figure
% Scaling the axes with increasing groove 
if(y4>=b_n_o)
    ylim(opt.axes_Animate_Nut,[-0.5,y4+0.5]);
    xlim(opt.axes_Animate_Nut,[-y4/2-0.5,y4/2+0.5]);
else
    ylim(opt.axes_Animate_Nut,[-0.5,x3*2+0.5]);
    xlim(opt.axes_Animate_Nut,[-x3-0.5,x3+0.5]);
end
ylabel(opt.axes_Animate_Nut,'[mm]')
xlabel(opt.axes_Animate_Nut,'[mm]')
title(opt.axes_Animate_Nut,'Nutmaﬂe')
hold(opt.axes_Animate_Nut, 'off')
