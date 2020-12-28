function S = nonLinearObs1(t,f,F,L)
M = 1000; 
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.8;

x = f(1);
dx = f(2);
t1 = f(3);
dt1 = f(4);
t2 = f(5);
dt2 = f(6);
out = [0; x; 0];
c1 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0];
control = L*(out - c1*f);
S=zeros(6,1);
S(1) = dx + control(1);
S(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (l1*m1*(S(3)^2)*sin(t1)) - (l2*m2*(S(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2))) + control(2);
S(3) = dt1 + control(3);
S(4) = ((cos(t1)*S(2)-g*sin(t1))/l1) + control(4);
S(5) = dt2 + control(5);
S(6) = (cos(t2)*S(2)-g*sin(t2))/l2 + control(6);
end