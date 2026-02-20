function [t,Y,num_fval] = Adapt_ERK_pq_STANDARD(f,tspan,y0,TOL,h0,butcher)
% Adaptive step size integrator using an embedded ERK method.

% INIZIALIZATION
t(1) = tspan(1);
t_curr = t(1);
h_curr = h0;
y_curr = y0;
Y(:,1) = y0; 

num_fval = 0;
p = butcher.p;
q = p+1;
k = 2;

% ADVANCING TIME WITH ADAPTIVE STEP
while t_curr < tspan(2)
    T_est = TOL + 1;
    cond=1;
    while cond
        
        %compute y_new with and T_est
        [y_new, ~, T_est] = RKq_Emb(f, t_curr, y_curr,  butcher,h_curr);
        num_fval = num_fval + p + 1;
        
        %acceptance: if the error is less then the tollarance
        if  T_est <= TOL  
            y_curr = y_new; 
            t_curr = t(k-1) + h_curr;
            cond=0;
            
        %in case of rejection, new step size
        else
            h_new = h_curr * abs(TOL / T_est)^(1 / (p+1));
            h_curr = h_new;
        end

    end
    
    % if %current time exceeds time's boundary
    if t_curr > tspan(2)
        t_curr = tspan(2);
        h_curr = t_curr - t(k-1);
        [y_new, ~,~] = RKq_Emb(f, t(k-1), Y(:,k-1), butcher, h_curr);
    end

    % store the solution and time
    Y(:,k) = y_new; 
    t(k) = t_curr;
    
    k = k+1;
end
end