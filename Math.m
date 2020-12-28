clc
clear

syms O_1 O_2 L1 L2 M m1 m2 g F;
syms x_dot O1_dot O2_dot

c1 = cos(O_1); s1 = sin(O_1); c2 = cos(O_2); s2 = sin(O_2);

%Velocity of mass 1 and 2 squared
v1 = (x_dot-L1*c1*O1_dot)^2 + (L1*s1*O1_dot)^2;
v2 = (x_dot-L2*c2*O2_dot)^2 + (L2*s2*O2_dot)^2;

%Kinetic Energy (T)
T = simplify(1/2*x_dot^2*M + 1/2*v1*m1 + 1/2*v2*m2);

%Potential Energy (V)
V = -m1*g*L1*c1 - m2*g*L2*c2;


%Lagrangian
L = T - V;
pretty(collect(L));



%Euler-Lagrange Equation
% dL/dX Term
dL__dx = 0;
dL__dO1 = diff(L,O_1);
dL__dO2 = diff(L,O_2);

syms O1(t) O2(t) x(t)

dL__dO1 = subs(dL__dO1, [x_dot O1_dot O2_dot O_1 O_2], [diff(x(t)) diff(O1(t)) diff(O2(t)) O1(t) O2(t)]);
dL__dO2 = subs(dL__dO2, [x_dot O1_dot O2_dot O_1 O_2], [diff(x(t)) diff(O1(t)) diff(O2(t)) O1(t) O2(t)]);


% d/dt(dL/dX_dot) Term
dL__dx_dot = diff(L,x_dot);
dL__dO1_dot = diff(L,O1_dot);
dL__dO2_dot = diff(L,O2_dot);

dL__dx_dot = subs(dL__dx_dot, [x_dot O1_dot O2_dot O_1 O_2], [diff(x(t)) diff(O1(t)) diff(O2(t)) O1(t) O2(t)]);
dL__dO1_dot = subs(dL__dO1_dot, [x_dot O1_dot O2_dot O_1 O_2], [diff(x(t)) diff(O1(t)) diff(O2(t)) O1(t) O2(t)]);
dL__dO2_dot = subs(dL__dO2_dot, [x_dot O1_dot O2_dot O_1 O_2], [diff(x(t)) diff(O1(t)) diff(O2(t)) O1(t) O2(t)]);

d_dt__dL__dx_dot = diff(dL__dx_dot,t);
d_dt__dL__dO1_dot = diff(dL__dO1_dot,t);
d_dt__dL__dO2_dot = diff(dL__dO2_dot,t);


% d/dt(dL/dX_dot) - dL/dX
Eq1 = simplify(d_dt__dL__dx_dot - dL__dx - F);
Eq2 = simplify(d_dt__dL__dO1_dot - dL__dO1);
Eq3 = simplify(d_dt__dL__dO2_dot - dL__dO2);


% Finding Equations of Motion
syms x_dot_dot O1_dot_dot O2_dot_dot

Eq1 = subs(Eq1,[diff(x(t),t,t) diff(O1(t),t,t) diff(O2(t),t,t)], [x_dot_dot O1_dot_dot O2_dot_dot]);
Eq2 = subs(Eq2,[diff(x(t),t,t) diff(O1(t),t,t) diff(O2(t),t,t)], [x_dot_dot O1_dot_dot O2_dot_dot]);
Eq3 = subs(Eq3,[diff(x(t),t,t) diff(O1(t),t,t) diff(O2(t),t,t)], [x_dot_dot O1_dot_dot O2_dot_dot]);

%Solving for intermediate theta_dot_dot equations
Theta1_dot_dot = solve(Eq2,O1_dot_dot);
Theta2_dot_dot = solve(Eq3,O2_dot_dot);

%Solving for Final x_dot_dot equation
Eq1 = simplify(subs(Eq1, [O1_dot_dot O2_dot_dot], [Theta1_dot_dot Theta2_dot_dot]));
Position_dot_dot = simplify(solve(Eq1,x_dot_dot));

Theta1_dot_dot = subs(Theta1_dot_dot, [x_dot_dot], [Position_dot_dot]);
Theta2_dot_dot = subs(Theta2_dot_dot, [x_dot_dot], [Position_dot_dot]);

pretty(simplify(Position_dot_dot))
pretty(simplify(Theta1_dot_dot))
pretty(simplify(Theta2_dot_dot))






