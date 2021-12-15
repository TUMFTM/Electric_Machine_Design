%*************************************************************************
% Problem Description:
%
% Reference : [1] Deb K, Pratap A, Agarwal S, et al. A fast and elitist 
%   multiobjective genetic algorithm NSGA-II[J]. Evolutionary Computation. 
%   2002, 6(2): 182-197.
% [2] T. Goel and N. Stander, “A Study on the Convergence of Multiobjective
% Evolutionary Algorithms,” Livermore Software Technology Corporation, 
% Hoboken, N.J, The Frank J. Fabozzi series, Nov. 2016.

%*************************************************************************

delete(gcp('nocreate'))

options = nsgaopt();                    % create default options structure
options.popsize = 80;                   % population/generation size, sollte Vielfaches der Anzahl der Kerne sein, Laptop hat 4, Cluster hat 32, default:40
options.maxGen  = 40;                  % max number of generation, wie viele iterationen, min: 10, default: 100

%% Add folders etc.
%addpath(genpath('W:\Projekte\Projekt DeTailED\intern\Modelle\24. MEAPA_Tool_v12 - Rechenknecht 2 - Kopie')); % anpassen!

different_color = true;

% Auswahl mit oder ohne GUI
Enable_GUI = 0;


%% Stopping Criteria [2] implemented in nsga2.m and ndsort.m

options.impr_ratio = 0.6;                % Minimum Improvement Ratio from gen to next gen
options.stop_crit = 'off';               % On oder Off

%% Set options for problem
options.numObj = 3;                     % number of objectives, y-Werte, angepasst auf 3 Werte
options.numVar = 10;                     % number of design variables, x-Werte
options.numCons = 0;                    % number of constraints
                  %P_n  n_N     U_n  p   Typ  Kue Magn  Batt Gear AnzMasch
options.vartype = [2    2       2    2    2   2   2     2     2     2];                  % 1: continous, 2: discrete
options.lb =      [70   7000    300  4    1   1   1     20    2     1];                  % lower bound of x
options.ub =      [200  10000  800  6    2   2   4     250   10     2];                  % upper bound of x
options.nameObj = {'Kosten in €'; 'Zykluseffizienz in %'; 'Maschinenvolumen'}; % Objectives benennen (Achsenbschriftung), erweitert auf 3 Werte

options.objfun = @optimization_Kalt_PMSM_objfun;     % objective function handle

options.plotInterval = 1;               % Song2011, S. 11 interval between two calls of "plotnsga". --> plot wird nach jeder iteration aufgerufen
options.outputInterval = 1;             % Song2011, S. 9, Intervall, welche Generationen er ablegt/speichert

%% Initialize first generation manually; Song2011, S. 6-7

%options.initfun={@initpop, strFileName, ngen}      % Restart from exist population file
%options.initfun={@initpop, oldresult, ngen}        % Restart from exist optimization result


%% Cossover-Options; Song2011, S. 7

%options.crossoverFraction = 0-1                     % Determines the share of the next generation produced by crossover. The remaining individuals of the next generation are produced by mutation. The crossover percentage must be between 0 and 1.
%options.crossover={'intermediate', ratio}           % ratio = “Verhältnis / Anteil”: child1 = parent1+rand*Ratio*(parent2 - parent1) 


%% Mutation - Options; Song2011, S. 8

%options.mutaionFraction = 1 - options.crossoverFraction     % Determines the share of the next generation produced by mutation. The remaining individuals of the next generation are produced by mutation. The mutation percentage must be between 0 and 1.
%options.mutation = {'gaussian', scale, shrink}              % scale: Deviation in first generation
                                                             % shrink: controls shrinkage of standard deviation as generations progress
%% Run parallel; Song2011, S. 11
% options.useParallel = 'no';             % parallel computation {‘yes’, ‘no’}
options.useParallel = 'yes';             % parallel computation {‘yes’, ‘no’}, mehrere Kerne benutzen
% options.poolsize = 32;                 % number of worker processes for FTM cluster
options.poolsize = 4;                    % number of worker processes for workstation, wieviele Kerne hat er

result = nsga2(options);                 % begin the optimization!

PEDS_Auswertung_Design_Variablen_Wolff;

%Animation_NSGAII_2D;
%Animation_NSGAII_3D;
