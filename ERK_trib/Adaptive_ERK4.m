function [t,u,h,num_fval] = Adaptive_ERK4(f,tspan,y0,h0,tol)

u(:,1)=y0; h(1)=h0; t(1)=tspan(1); T=tspan(2); num_fval=0;

% Bogacki-Shampine formula
c=[0,1/2,3/4,1]; 
A=[0 0 0 
   1/2 0 0
   0 3/4 0 
   2/9 1/3 4/9];
bp=[2/9 1/3 4/9 0];   
bq=[7/24 1/4 1/3 1/8];

k=1; t_curr=t(1); T_est=tol+1;
while t_curr<T
    k=k+1;
    h_curr=h(k-1);
    
    while  T_est >= tol
        num_fval=num_fval+4;
        K1=f(      t(k-1),               u(:,k-1) );
        K2=f(   t(k-1)+c(2)*h_curr,       u(:,k-1)+h_curr*A(2,1)*K1 );
        K3=f(   t(k-1)+c(3)*h_curr,   u(:,k-1)+h_curr*A(3,1)*K1+h_curr*A(3,2)*K2 );
        K4=f(   t(k-1)+c(4)*h_curr,   u(:,k-1)+h_curr*A(4,1)*K1+h_curr*A(4,2)*K2+h_curr*A(4,3)*K3);     
        sum_p=(bp(1)*K1+bp(2)*K2+bp(3)*K3);
        T_est=  h_curr * norm(  sum_p - (bq(1)*K1+bq(2)*K2+bq(3)*K3+bq(4)*K4) );
        if T_est >= tol
            h_curr=h_curr* norm(tol/T_est)^(1/3);
        end
    end
    K1=f(      t(k-1),               u(:,k-1) );
    K2=f(   t(k-1)+c(2)*h_curr,       u(:,k-1)+h_curr*A(2,1)*K1 );
    K3=f(   t(k-1)+c(3)*h_curr,   u(:,k-1)+h_curr*A(3,1)*K1+h_curr*A(3,2)*K2 );     
    num_fval=num_fval+3;
    sum_p=(bp(1)*K1+bp(2)*K2+bp(3)*K3);
    h(k,1)=h_curr;
    u(:,k)=u(:,k-1)+h(k)*sum_p;
    t(k,1)=t(k-1)+h(k); 
    t_curr=t(k);
end
end