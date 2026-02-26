syms t y(t)
f=@(t,y)(3+cos(t)-y);
ode= diff(y,t)==3+cos(t)-y; cond= y(0)==1;
ysolv(t)=dsolve(ode,cond);
sol1=matlabFunction(ysolv);
figure(99)
fplot(ysolv,[0 20])
title("sol of y'= f(t,y)")

% dati relativi al pb
N=100;  y0=1; tspan=[0 20]; 
h=(tspan(2)-tspan(1))/N;


% ERK3
b1=1/4; b2=0;
[t,u1,num_fval,b,c,A] = ERK3_2inputs(f,tspan,y0,h,b1,b2);
figure(1)
subplot(2,2,2)
plot(t,u1)
title("ERK3 with b1=1/4, b2=0")

b1=1/6; b2=2/3;
[t,u2,num_fval,b,c,A] = ERK3_2inputs(f,tspan,y0,h,b1,b2);
figure(1)
subplot(2,2,3)
plot(t,u2)
title("ERK3 with b1=1/6, b2=2/3")

b1=1/4; b2=1/4;
[t,u3,num_fval,b,c,A] = ERK3_2inputs(f,tspan,y0,h,b1,b2);
figure(1)
subplot(2,2,4)
plot(t,u3)
title("ERK3 with b1=1/4, b2=1/4")

b=[1/4,0,3/4]; c=[0,1/3,2/3]; 
A=[0 0
   1/3 0
   0 2/3];
[t,u4,num_fval] = ERK3(f,tspan,y0,h,b,c,A);
figure(1)
subplot(2,2,1)
plot(t,u4)
title("ERK3 with Heun coeff.")
sgtitle("3-stage Runge-Kutta method")

% ERK4
b=[1/8,3/8,3/8,1/8]; c=[0,1/3,2/3,1]; 
A=[0 0 0
   1/3 0 0 
   -1/3 1 0 
   1 -1 1];
[t,u5,num_fval] = ERK4(f,tspan,y0,h,b,c,A);
figure(2)
subplot(1,2,1)
plot(t,u5)
title("ERK4 with first set of coeff.")

b=[1/6,1/3,1/3,1/6]; c=[0,1/2,1/2,1]; 
A=[0 0 0
   1/2 0 0 
   0 1/2 0 
   0 0 1];
[t,u6,num_fval] = ERK4(f,tspan,y0,h,b,c,A);
figure(2)
subplot(1,2,2)
plot(t,u6)
title("ERK4 with second sett of coeff.")
sgtitle("4-stage Runge-Kutta method")

% Adaptive ERK4
h0=0.3; tol=10^-2;
[t,u7,h,num_fval] = Adaptive_ERK4(f,tspan,y0,h0,tol);
figure
plot(t,u7)
title("Adaptive ERK4")





