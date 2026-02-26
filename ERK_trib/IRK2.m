function [t,u,num_fval] = IRK2(f,tspan,y0,h,butcher,tau,maxit)
b=butcher.b; A=butcher.A; c=butcher.c;
q = length(b);  % times of steps to do
n = length(y0);
u(:,1)=y0; 
t=(tspan(1):h:tspan(2));
N=length(t); 
num_fval=0; 
K0=ones(q*n,1);
for i=2:N
    tcurr=t(i-1); 
    ucurr=u(:,i-1);
    F = @(x) [
        x(1:n) - f(tcurr + h * c(1), ucurr + h * (A(1,1) * x(1:n) + A(1,2) * x(n+1:2*n)));
        x(n+1:2*n) - f(tcurr + h * c(2), ucurr + h * (A(2,1) * x(1:n) + A(2,2) * x(n+1:2*n)))
    ];
    [K, N_fval] = NWT_solver(K0, F, tau, maxit);
    K=reshape(K,[n,q]);
    Z=zeros(n,1);
    for j=1:q
        Z=Z+b(j)*K(:,j); 
    end
    u(:,i)=u(:,i-1)+h*Z;
    num_fval=num_fval+N_fval;   
end
end


   