function [t,Y,num_fval] = Fixed_ERK_q(f,tspan,y0,h0,butcher)

t(1) = tspan(1);
Y(:,1) = y0;

num_fval = 0;

A = butcher.A; 
b = butcher.b; 
c = butcher.c;
q = length(butcher.c);

K=[];

N = round((tspan(2)-tspan(1))/h0); 

%COMPUTE THE SOLUTION FOR EACH TIME STEP
for i=2:N+1

    % COMPUTE 1 STAGE
    K(:,1) = f(t(i-1), Y(:,i-1));

    %COMPUTE FOLLOWING STAGES
    for j=2:q
        K(:,j) = f(t(i-1) + c(j)*h0, Y(:,i-1) + h0*(K(:,1:j-1)*A(j,1:j-1)')); 
    end

    %STORE THE SOLUTION
    Y(:,i) = Y(:,i-1) + h0*(K*b);
    t(i) = t(i-1) + h0;
    
    num_fval = num_fval + q;

end