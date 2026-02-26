N=200; 
tau=10^-3; maxit=100;
c=[(3-sqrt(3))/6 (3+sqrt(3))/6]';
b=[1/2 1/2]';
A=[1/4 (3-2*sqrt(3))/12
  (3+2*sqrt(3))/12 1/4];
butcher=struct("A",A,"b",b,"c",c);

% funzione a valori in R
f=@(t,y)(3+cos(t)-y); y0=1; tspan=[0 20]; 
h=(tspan(2)-tspan(1))/N;

[t,u,~] = IRK(f,tspan,y0,N,butcher);
IRK2(f,tspan,y0,h,butcher,tau,maxit)
figure(99)
subplot(1,2,1)
plot(t,u)
title("irk")
[t,u,~] = IRK2(f,tspan,y0,h,butcher,tau,maxit);
figure(99)
subplot(1,2,2)
plot(t,u)
title("irk 2 steps")

% funzione a valori in Rn
y0=ones(3,1); tspan=[0 100];
Ti=linspace(tspan(1),tspan(2),N); h=(tspan(2)-tspan(1))/N;
k1=1; k2=10;
f=@(t,y)[-k1*y(1);
          k1*y(1)-k2*y(2);
          k2*y(2)];
[t,u,h] = IRK(f,tspan,y0,N,butcher);
figure(1)
plot(t,u)
title("irk")
