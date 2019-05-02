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

function [i_dsol, i_qsol] = Optimierung_i(prim, ent, reg, omega_mesh, M_mesh)
%Optimization minimizes the d- and q-component of the current under
% secondary conditions and limits
% Notation: i_d = i(1), i_q = i(2)
% Bound: -inf<=i_d<=0, 0<=i_q<=inf
% Linear Inequality Constraint: none
% Linear Equality Constraint: none
% Nonlinear constraints: see function nonlcon_fun
% For reasons of calculation time, values from column 2 or row 2 onward
% of previous currents are used as start values.

%% Preprocessing
% Minimization function 
fun = @(i) i(1)^2 + i(2)^2;
% initial value
x0 = [0 0];
% Bound Constraints
lb = [-inf 0];
ub = [0 inf];
% Linear Inequality Constraint
A = [];
b = [];
% Linear Equality Constraint
Aeq = [];
beq = [];
% Nonlinear Constraints
nonlcon = @nonlcon_fun;
% Options for Optimization
options = optimoptions(@fmincon,'Display','off','Algorithm','sqp');

% Allocate memory
i_dsol = zeros(size(M_mesh));
i_qsol = zeros(size(M_mesh));

%% Solving
tic
for j = 1:length(omega_mesh(1,:))
    omega_k = omega_mesh(1,j);
    for k = 1:length(M_mesh(:,1))
        M = M_mesh(k,j);
        if(isnan(M))
            % Sorting out points that cannot be reached (full load characteristic)
            i_dsol(k,j) = NaN;
            i_qsol(k,j) = NaN;
        else
            if(j>1 && ~isnan(i_dsol(k,j-1)) && ~isnan(i_qsol(k,j-1)))
                if(k>1 && ~isnan(i_dsol(k-1,j)))
                    % Minimization with initial values from previous calculation
                    [x,~,exitflag,~] = fmincon(fun, [i_dsol(k-1,j) i_qsol(k,j-1)], A, b, Aeq, beq, lb, ub, ...
                    @(i) nonlcon(i, M, reg.u_max, reg.i_max, ent.L_d, ent.L_q, ent.R_s, omega_k, ent.psi_PM, prim.p), options);
                else
                    % Minimization with initial values from previous calculation
                    [x,~,exitflag,~] = fmincon(fun, [i_dsol(k,j-1) i_qsol(k,j-1)], A, b, Aeq, beq, lb, ub, ...
                        @(i) nonlcon(i, M, reg.u_max, reg.i_max, ent.L_d, ent.L_q, ent.R_s, omega_k, ent.psi_PM, prim.p), options);
                end
            else
                % Minimization with initial value [0,0]
                [x,~,exitflag,~] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, ...
                @(i) nonlcon(i, M, reg.u_max, reg.i_max, ent.L_d, ent.L_q, ent.R_s, omega_k, ent.psi_PM, prim.p), options);
            end
            % Values < tol are set to 0 (numerical inaccuracy)
            if(abs(x(1))<1e-2)
                x(1)=0;
            end
            if(abs(x(2))<1e-2)
                x(2)=0;
            end
            if(exitflag~=1 && exitflag~=2)
                % Current is set to inf if no result is found
                i_dsol(k,j) = inf;
                i_qsol(k,j) = inf;
            else
                % Assign values when result is found
                i_dsol(k,j) = x(1);
                i_qsol(k,j) = x(2);
            end
        end
    end
end
toc

end

function [c,ceq] = nonlcon_fun(i, M_ref, u_max, i_max, L_d, L_q, R_s, omega_k, psi_PM, p)
% nonlcon_fun Definition of the nonlinear constraints
% Secondary condition 1: Current limit
% Secondary condition 2: Voltage limit
% Secondary condition 3: Moment reference
    
    c = [i(1)^2 + i(2)^2 - i_max^2;...
        (R_s*i(1)-omega_k*L_q*i(2))^2 + (R_s*i(2)+omega_k*L_d*i(1)+omega_k*psi_PM)^2 - u_max^2];
    ceq = 1.5*p*(L_d-L_q)*i(2)*i(1) + 1.5*p*psi_PM*i(2) - M_ref;
end