clc

syms M m1 m2 g l1 l2 x xd o1 o1d o1dd o2 o2d o2dd F
syms f1 f2 f3 f4 f5 f6 

xdd_n = F- m1*g*sin(o1)*cos(o1)- m2*g*sin(o2)*cos(o2)- m1*l1*(o1d^2)*sin(o1)- m2*l2*(o2d^2)*sin(o2);
xdd_d = (M+m1+m2) - m1*(cos(o1))^2 - m2*(cos(o2))^2;
xdd = xdd_n/xdd_d;

o1dd= (xdd*cos(o1)-g*sin(o1))/l1;

o2dd= (xdd*cos(o2)-g*sin(o2))/l2;

f1= xd;
f2= xdd;
f3= o1d;
f4= o1dd;
f5= o2d;
f6= o2dd;

f = [f1 f2 f3 f4 f5 f6];
q = [x xd o1 o1d o2 o2d];



A = cell(6); %cell(6,6);
x= 0;
o1= 0;
o2= 0;
o1d = 0;
o2d = 0;
xd = 0; 



B = cell(1, 6);
for i=1:6
    f_ = f(i);
    B{i} = diff(f_, F);
    for j=1:6
        q_ = q(j);
        diff_ans = diff(f_, q_);
               
        subs_ans = subs(diff_ans);
        A{i,j} = subs_ans;      
    end
end

M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 10;

% Jacobian - A and B
A_subs = zeros(6);
B_subs = zeros(1,6);
for i=1:6
    B_subs(i) = subs(B(i));
    for j=1:6
        A_subs(i,j) = subs(A{i,j});      
    end
end


% controbility check
CC = ctrb(A_subs, B_subs');

rank(CC);

Q = zeros(6);

Q(1,1) = 10;
Q(2,2) = 1;
Q(3,3) = 150000;
Q(4,4) = 1;
Q(5,5) = 150000;
Q(6,6) = 1;

R = 1;

poles = eigs(A_subs);
% getting value of gain by lqr
K = lqr(A_subs,B_subs',Q,R) 
% % setting initial position
% X0 = [5 0 pi/6 0 pi/12 0]
 D = [0];
% % implementing lqr
Ac = (A_subs - B_subs'*K);
poles_c = eigs(Ac)
% C1 = [1 0 0 0 0 0; 0 0 1 0 0 0;0 0 0 0 1 0];
% lqr_sys = ss(Ac,B_subs',C1, D);
% time= 0:0.1:500
% [y,t,x]=initial(lqr_sys,X0,time);
% figure,
% plot(t, y(:,1));
% xlabel('Time');
% ylabel('Position');
% 
% figure,
% plot(t, y(:,2));
% xlabel('Time');
% ylabel('Pendulum1 angle');
% 
% figure,
% plot(t, y(:,3));
% xlabel('Time');
% ylabel('Pendulum2 angle');

%plot(t,y);
C1 = [1 0 0 0 0 0];
C2 = [0 0 1 0 0 0; 0 0 0 0 1 0];
C3 = [1 0 0 0 0 0; 0 0 0 0 1 0];
C4 = [1 0 0 0 0 0; 0 0 1 0 0 0;0 0 0 0 1 0];
Ob_po = obsv(A_subs,C1)
disp('rank is');
rank(Ob_po)

% place_poles = [2i , -12, -2i, -7, -3, -4];
 place_poles = [-11 , -12, -10, -7, -3, -4];

%implmenting luenberger

% L1 = place(A_subs', C1', place_poles)' 
% luenburger1 = ss(A_subs - L1*C1, B_subs', C1, D)
% step(luenburger1)

% L2 = place(A_subs', C2', place_poles)' 
% luenburger2 = ss(A_subs - L2*C2, B_subs', C2, D)
% step(luenburger2)
% 
% L3 = place(A_subs', C3', place_poles)' 
% luenburger3 = ss(A_subs - L3*C3, B_subs', C3, D)
% step(luenburger3)
% 
% 
L4 = place(A_subs', C4', place_poles)';
X0 = [5 0 pi/6 0 pi/12 0];
luenburger4 = ss(A_subs - L4*C4, B_subs', C4, D);
step(luenburger4);

% % Kalman Estimator 
% tspan = 0:0.1:100;
% B = 0.1*eye(6); %Process Noise
% V = 0.01; %Measurement Noise
% c = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0];
% L = lqe(A,B,c,B,V*eye(3));
% Ac1 = A-(L*c);
% kalman_sys = ss(Ac1,[B L],c,0);
% 
% % Non-linear LQG Response
% q_initial = [5 0 pi/6 0 pi/12 0];
% [t,q1] = ode45(@(t,q)NonLinear_LQG(t,q,-K*q, L),tspan,q_initial);
% figure();
% hold on
% plot(t,q1(:,5))
% ylabel('q_state')
% xlabel('time')
% legend('x')
% hold off



