clear; clc;
syms x1 x2;

% declare the state
vars = [x1; x2];

% declare the system dynamics
f = [-x1 + x2^3;
     -x2];

% create SOSTOOLS program
prog = sosprogram(vars);

% Lyapunov is defined second-order
monom = monomials(vars, 0:2);  % monomial basis e.g., 2 --> x1^2, x1*x2, x2^2
[prog, V] = sospolyvar(prog, monom, 'wscoeff');

% V(x) > ε‖x‖²
eps = 1e-4;
prog = sosineq(prog, V - eps*(x1^2 + x2^2));  % guarantee V > 0 and V([0,0]) == 0

% calculate the derivative of V
Vdot = jacobian(V, vars) * f;

% negative definite guarantee: -Vdot > 0, Vdot < 0
prog = sosineq(prog, -Vdot);

% solve
prog = sossolve(prog);

% result
Vsol = sosgetsol(prog, V);
disp('Lyapunov function V(x):')
disp(Vsol)
