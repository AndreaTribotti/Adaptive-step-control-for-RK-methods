function [t,u,num_fval] = ERK3(f,tspan,y0,h,b,c,A)
u(:,1)=y0; t(1)=tspan(1);
N=round((tspan(2)-tspan(1))/h);
for k=2:N+1
    K1=f(      t(k-1),               u(:,k-1) );
    K2=f(   t(k-1)+c(2)*h,       u(:,k-1)+h*A(2,1)*K1 );
    K3=f(   t(k-1)+c(3)*h,   u(:,k-1)+h*A(3,1)*K1+h*A(3,2)*K2 );
    u(:,k)=u(:,k-1)+h*(b(1)*K1+b(2)*K2+b(3)*K3);
    t(k,1)=t(k-1)+h;
end
num_fval=3*N;
end