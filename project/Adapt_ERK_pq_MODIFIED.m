function [t,Y,num_fun_eval, H_acp, H_rej, t_rej, local_err] = Adapt_ERK_pq_MODIFIED(f,tspan,y0,flag,butcher,Atol,Rtol)
% embedded Runge-Kutta with adaptive step size

% INITIAL CONDITIONS
Y(:,1) = y0;                    
y_curr = y0;                    
t(1) = tspan(1);                
t_curr = t(1);         

num_fun_eval = 2;              
p = butcher.p; % principal order
q = p+1;         
k = 2;

f_curr = f(t_curr, y_curr); 
n = length(f_curr);

H_acp = []; % accepted step sizes
H_rej = []; % rejected step sizes
t_rej = []; % moments of the rejection

local_err = [];  %local errors

%Safety factor
if flag==0
    fac = 0.8;   
elseif flag ==1
    fac = 0.9;
elseif flag ==2
    fac = (0.25)^(1/p+1);   
elseif flag == 3
    fac = (0.38)^(1/p+1); 
end

% Safety factor's span
facmax = 3.8;         
facmin = 0.5;        

%Max number of iterations:
max_iter = 80;

% STARTING STEP SIZE
% previous and embedded solution are zero vectors
d0 = err_norm(zeros(n,1),y0,zeros(n,1),Atol,Rtol);           
d1 = err_norm(zeros(n,1),f_curr,zeros(n,1),Atol, Rtol);     

if d0 < 1e-5 || d1 < 1e-5
    h0 = 1e-6;      
else
    h0 = 0.01 * d0 / d1;
end

y1 = y0 + h0*f_curr;
f1 = f(t_curr + h0, y1);

d2 = (err_norm(zeros(n,1),f1,f_curr,Atol, Rtol))/h0;

if max(d1,d2) <= 1e-15
    h1 = max(1e-6, h0*1e-3);
else
    h1 = (0.01 \ max(d1,d2)) ^ (1 / (p+1));
end

h_curr = min(100*h0, h1);


% ADAPTATION
while t_curr < tspan(2)
    cond = 1;       
    num_iter = 0;

    while cond 
        % computing y_new with h_old using the embedded RK step
        [y_new,y_emb,T_est] = RKq_Emb(f, t_curr, y_curr, butcher, h_curr);
        [err] = err_norm(y_curr,y_new,y_emb,Atol,Rtol);
        num_fun_eval = num_fun_eval + q; 
        num_iter = num_iter + 1;

        %acceptance criteria
        if  err<=1  && (num_iter > max_iter  || 0.5<err)            
            y_curr = y_new; 
            t_curr = t(k-1) + h_curr;
            cond = 0;

        %in case of rejection, new step size
        else
            h_new = h_curr * min(facmax, max(facmin, fac * abs(1/err).^(1 /(p+1))));
            h_curr = h_new;

            H_rej(end+1)=h_curr;
            t_rej(end+1)=t_curr;
        end

    end

    % if current time exceeds time's boundary
    if t_curr > tspan(2)
        t_curr = tspan(2);
        h_curr = t_curr - t(k-1);
        [y_new, ~, T_est] = RKq_Emb(f, t(k-1), Y(:,k-1), butcher, h_curr);
        num_fun_eval = num_fun_eval + q;
    end

    % store the accepted solution and time
    Y(:,k) = y_new; 
    t(k) = t_curr;
   
    k = k+1;

    %estimated local error
    local_err(end+1)=T_est;


end

H_acp=t(2:end) - t(1:end-1);

end