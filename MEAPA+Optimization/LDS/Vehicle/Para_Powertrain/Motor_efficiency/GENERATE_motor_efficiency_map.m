function [motor] = GENERATE_motor_efficiency_map(motor)

% This function generates a motor efficiency map of an ASM or PMSM with 
% respect to the specified top level design parameters 

%% Initialize from input
% Defined top level parameters
motortype = motor.type;     % PMSM or ASM
P_n = motor.power;          % in kW
n_n = motor.n_n;         	% in 1/min
n_max = motor.n_max;        % in 1/min
U_n = motor.U_n;            % in V

% Specific parameter for ASM (fixed!)
p = 2;                   % number of pole pairs
m= 3;                    % number of phases

% Specific parameter for PMSM (fixed!)
Uberlastfaktor = 1.0;      % overload factor
    
%% path names
function_name = '\GENERATE_motor_efficiency_map';
base_path = mfilename('fullpath');
base_path = erase(base_path,function_name);
warning ('off','all');

%% Number of data points of the efficiency map
nn = 160; %number of speed data points 
nT = 160; %number of torque data points 

%% Check if new motor map is required
mot_path = [motor.type '_' num2str(P_n) '_' num2str(n_n) ...
    '_' num2str(n_max) '_' num2str(U_n) ]; %name of motor file

IsNewBuild = exist([base_path '\MOT_memory\' mot_path '.mat']) ~= 2;

if IsNewBuild == 0 %if motor already exists, just load
    load([base_path '\MOT_memory\' mot_path '.mat']);

else %if not, create new motor
%% Efficiency map generation
if strcmp(motortype, "ASM")
    %% Add and delete paths to avoid erronous usage of functions with same name 
    addpath(genpath([base_path '\ASM']));
    rmpath(genpath([base_path '\PMSM']));

    %% Calculation
    [eff,eff_n_axis,eff_T_axis,GK_M_W_Matrix_6] = ASM_Kennfeld_Tool(P_n,n_n,n_max,p,U_n,m);
    close all

    eff = eff(length(eff_T_axis):end,:); %take only half of the map (--> motor part)
    eff_n_axis = eff_n_axis*60/2/pi; % in 1/min

    T_max = GK_M_W_Matrix_6(5,:); % Maximum torque curve
    T_max_n_axis = eff_n_axis; % speed axis of maximum torque curve

    clear GK_M_W_Matrix_6  

elseif strcmp(motortype, "PMSM")
    %% Add and delete paths to avoid erronous usage of functions with same name 
    addpath(genpath([base_path '\PMSM']));
    rmpath(genpath([base_path '\ASM']));

    % Specific parameter for PMSM (fixed!)
    P_n = P_n*1e3;           % in W - nominal power

    %% Calculation
    [om_M_eta,OM_Vektor,M,m] = g00_Synchronmotorauslegung(P_n,n_n,n_max,U_n,Uberlastfaktor); 

    %% Results
    eff = om_M_eta;        
    eff_n_axis = OM_Vektor*60/2/pi;
    eff_T_axis = M(1,:);   

    T_max = m; % Maximum torque curve
    T_max_n_axis = eff_n_axis; % speed axis of maximum torque curve 

    %% Allocate to MOT struct to manipulate max. torque and speed if they are not correct
    MOT.eff = eff; 
    MOT.eff_n_axis = eff_n_axis;
    MOT.eff_T_axis = eff_T_axis;
    MOT.T_max = T_max;
    MOT.T_max_n_axis = T_max_n_axis;
    
    %if the T_max curve is incorret: derive from eff_map
    if m(1) == m(end); [MOT] = derive_T_max_mot(MOT); end
    % trim and scale eff map to correct max. speed
    MOT = trim_and_scale(MOT, n_max);
    
    %Allocate back to have the same after processing like the ASM
    eff = MOT.eff; 
    eff_n_axis = MOT.eff_n_axis;
    eff_T_axis = MOT.eff_T_axis;
    T_max = MOT.T_max;
    T_max_n_axis = MOT.T_max_n_axis;
    
    clear om_M_eta OM_Vektor M m
else
    disp("An unexpected motortype was selected. Please check the input")  
end

if IsNewBuild == 1
    %% Interpolate efficiency map
%     [eff_n_axis, eff_T_axis, eff] = ...
%         Regrid_map(eff_n_axis,eff_T_axis,eff,nn,nT,'linear'); %interpolieren

    % Manipulate efficiency map
    eff(eff<0.1) = 0.1; %set small efficiency values to 0.1
    eff(1,:) = 0; %first colomn is zero
    eff(:,1) = 0.1; %first row is zero
    % Create one "larger" eff. map that includes acceleration and recuperation
    eff_recu = 1./eff; %recuperation efficiency
    eff_recu(eff_recu==inf)=100;
    eff_recu = eff_recu(2:end,:);
    eff = [flipud(eff_recu); eff ];
    eff_T_axis = [fliplr(-1.*eff_T_axis(2:end)) eff_T_axis];
    
    % maximum torque curve
    F = griddedInterpolant(T_max_n_axis,T_max,'linear','none');
    T_max_n_axis = 0: max(T_max_n_axis)/(nn-1) :max(T_max_n_axis);
    T_max = F(T_max_n_axis);

%% Allocate the results
% Efficiency map
motor.eff_n_axis    = eff_n_axis;
motor.eff_T_axis    = eff_T_axis;
motor.eff           = eff;
motor.eff_recu      = eff_recu;
% Maximum torque curve
motor.T_max         = T_max;
motor.T_max_n_axis  = T_max_n_axis;

%% Visualize (to check if maps are correct)
% figure;
% contLines = [0.6 0.8 0.85 0.9:0.02:1];
% contLines = [contLines 1./contLines];
% [C,h] = contourf(motor.eff_n_axis,motor.eff_T_axis,motor.eff,contLines);
% clabel(C,h); colorbar; clear C h;
% 
% figure; 
% plot(motor.T_max_n_axis,motor.T_max);

%% Save motor to load instead of calculate a new one
if motor.type == 'PMSM'; P_n = P_n/1000; end %motor power in kW
mot_path = [motor.type '_' num2str(round(P_n,2)) '_' num2str(n_n) ...
    '_' num2str(n_max) '_' num2str(U_n) ]; %name of motor file

save([base_path '\MOT_memory\' mot_path ], 'motor');

end

end
% Turn on warnings again
warning ('on','all');

end


%% derive max torque curve from eff. profile
function [MOTOR] = derive_T_max_mot(MOTOR)
%Help function to derive the maximum torque curve from the efficiency map,
%because the Horlbeck tool delivers an incorrect max. mot. torque curve for
%some input parameters
    
    %search all colomns
    [r,c] = size(MOTOR.eff);
    T_max = zeros(1,r);
    for i = 1:r
        if sum(isnan(MOTOR.eff(i,:)))>=1
            [r_nan] = find(isnan(MOTOR.eff(i,:)));
            T_max(i) = min(r_nan)*max(MOTOR.T_max)/c;
        else
            T_max(i) = max(MOTOR.T_max);
        end
    end
    MOTOR.T_max = T_max;
end

function [MOTOR] = trim_and_scale(MOTOR, n_max_desired)
%The Horlbeck tool calculates motor maps that have too large n_max (even if
%defined differently). A comparison with real motor maps shows, that the
%data is reasonable for a wide range of the desired field of operation.
%Hence, the map is cut at the desired maximum speed and re-scaled

    n = MOTOR.eff_n_axis(MOTOR.eff_n_axis<=n_max_desired);
    T = MOTOR.eff_T_axis;
%   trim
    map = MOTOR.eff(1:length(n),:);
    T_max = MOTOR.T_max(1:length(n));
%	rescale eff. map
    [nx,ny] = size(MOTOR.eff);
    [X, Y, MAP] = Regrid_map(n,T,map',nx,ny,'none');
    MOTOR.eff_n_axis = X;
    MOTOR.eff_T_axis = Y;
    MOTOR.eff = MAP;
%   rescale max torque curve 
    F = griddedInterpolant(n,T_max);
    MOTOR.T_max = F(X);
    MOTOR.T_max_n_axis = X;
end
