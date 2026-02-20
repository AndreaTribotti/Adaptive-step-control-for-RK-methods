function [y_HO,y_EMB,T_est] = RKq_Emb(f, t_curr, y_curr, butcher,  h)

A = butcher.A; 
b = butcher.b; 
c = butcher.c; 
q = length(c);
K=[];

% COMPUTE 1 STAGE 
K(:,1) = f(t_curr, y_curr);

%COMPUTE REMAINING STAGES
for i=2:q
    K(:,i) = f(t_curr + c(i)*h, y_curr + h*(K(:,1:i-1)*A(i,1:i-1)'));
end

%HIGHER ORDER SOLUTION (using 1 row of b)
y_HO = y_curr + h*(K*b(:,1));

%LOWER ORDER SOLUTION-EMBEDDED (using 2 row of b)
y_EMB = y_curr + h*(K*b(:,2));

%ESTIMATED TOLLERANCE
T_est = norm((h * (K* ( b(:,1)-b(:,2) ) )));
end