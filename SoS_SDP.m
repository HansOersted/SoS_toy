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
monom = monomials(vars, 2);  % monomial basis（x1^2, x1*x2, x2^2）
[prog, V] = sospolyvar(prog, monom, 'wscoeff');

% V(x) > ε‖x‖²
eps = 1e-4;
prog = sosineq(prog, V - eps*(x1^2 + x2^2));  % guarantee V > 0

% calculate the derivative of V
Vdot = jacobian(V, vars) * f;

% negative definite guarantee：-Vdot > ε‖x‖², Vdot < -ε‖x‖²
prog = sosineq(prog, -Vdot - eps*(x1^2 + x2^2));

% solve
prog = sossolve(prog);

% result
Vsol = sosgetsol(prog, V);
disp('Lyapunov function V(x):')
disp(Vsol)
