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

function [i_s, beta] = Optimierung_M(prim, ent, reg, omega_mesh)
 %Optimization maximizes torque under constraints and limits
% Notation: i_s = x(1), beta = x(2)
% Bound: 0<=i_s<=i_max, 0<=beta<=pi
% Linear Inequality Constraint: none
% Linear Equality Constraint: none
% Nonlinear constraints: see function nonlcon_fun

%% Preprocessing
% Minimization function
fun = @(x) -1.5.*prim.p.*((ent.L_d-ent.L_q).*x(1).^2.*0.5.*sin(2.*x(2)) + ent.psi_PM.*x(1).*sin(x(2)));
% initial value
x0 = [0 0];
% Bound Constraints
lb = [0 0];
ub = [reg.i_max pi];
% Linear Inequality Constraint
A = [];
b = [];
% Linear Equality Constraint
Aeq = [];
beq = [];
% Nonlinear Constraints
nonlcon = @nonlcon_fun1;
% Options for Optimization
options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point');

% Allocate memory
i_s = zeros(length(omega_mesh(1,:)),1);
beta = zeros(length(omega_mesh(1,:)),1);

%% Solving
tic
for j = 1:length(omega_mesh(1,:))
    omega_k = omega_mesh(1,j);
    % Minimization of initial value [0,0]
    [x_sol,~,exitflag,~] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, ...
    @(x) nonlcon(x, reg.u_max, ent.L_d, ent.L_q, ent.R_s, omega_k, ent.psi_PM), options);
    if(exitflag~=1 && exitflag~=2)
        % Current is set to inf if no result is found
        i_s(j) = inf;
        beta(j) = inf;
    else
        % Assign values when result is found
        i_s(j) = x_sol(1);
        beta(j) = x_sol(2);
    end
end
toc

end

function [c,ceq] = nonlcon_fun1(x, u_max, L_d, L_q, R_s, omega_k, psi_PM)
%nonlcon_fun Definition of nonlinear constraints
% Secondary condition 1: Voltage limit
    
    c = (R_s*x(1)*cos(x(2)) - omega_k*L_q*x(1)*sin(x(2)))^2 + (R_s*x(1)*sin(x(2)) + omega_k*L_d*x(1)*cos(x(2)) + omega_k*psi_PM)^2 - u_max^2;
    ceq = [];
end