%% SCRIPT OF THE ADAPTIVE STEP SIZE PROJECT
% We will compare the implementation of the adaptive step proposed in the
% pdf with the one we saw in class and the fixed one to see how the results changes 

clear
close all
load adapt_butchers.mat
load ERK4_butchers.mat

% SET THE EQUATION: BRUSSELATOR PROBLEM
tspan = [0 20];
y0 = [1.5; 3];
f = @(x,y) ([1 + y(1)^2 * y(2) - 4*y(1); ...
             3 * y(1) - y(1)^2 * y(2)] );
Atol = 1e-4; 
Rtol = 1e-4;

h0 = 1;

%% SOLVE THE BRUSSELATOR PB WITH DIFFERENT METHODS

% 0. FIXED SIZE STEP METHOD
[t_FX,Y_FX,num_fun_eval_FX] = Fixed_ERK_q(f,tspan,y0,h0,butcher1_ERK4);



% 1. STANDARD ADAPTIVE METHOD
[t_ST,Y_ST,num_fun_eval_ST] = Adapt_ERK_pq_STANDARD(f,tspan,y0,Rtol,h0,adapt_Kutta_Fehlberg_4_5);


% 2. MODIFIED ADAPTIVE METHOD VARYING THE FLAG PARAMETER:
%   --> changing the safety factor 

[t_M0,Y_M0,num_fun_eval_M0,H0_acp, H0_rej, t0_rej, local_err_0] = Adapt_ERK_pq_MODIFIED(f,tspan,y0,0,adapt_Kutta_Fehlberg_4_5,Atol,Rtol);
[t_M1,Y_M1,num_fun_eval_M1,H1_acp, H1_rej, t1_rej, local_err_1] = Adapt_ERK_pq_MODIFIED(f,tspan,y0,1,adapt_Kutta_Fehlberg_4_5,Atol,Rtol);
[t_M2,Y_M2,num_fun_eval_M2,H2_acp, H2_rej, t2_rej, local_err_2] = Adapt_ERK_pq_MODIFIED(f,tspan,y0,2,adapt_Kutta_Fehlberg_4_5,Atol,Rtol);
[t_M3,Y_M3,num_fun_eval_M3,H3_acp, H3_rej, t3_rej, local_err_3] = Adapt_ERK_pq_MODIFIED(f,tspan,y0,3,adapt_Kutta_Fehlberg_4_5,Atol,Rtol);



%% COMPARE THE MODELS:
% We will compare the models writing the relative error of both wrt the
% standard ode solver

%We will store  in a vector the datas characterizing this model:
%     Num of time steps
%     Max and min step sizes
%     Relative error
%     Num of function evaluations.


% 0. FIXED STEP SIZE METHOD
[~,Y_FX_True] = ode15s(f, t_FX, y0);
    
rel_err_FX = norm(Y_FX_True' - Y_FX) / norm(Y_FX_True');
Data_FX = [length(t_FX), h0, h0, rel_err_FX, num_fun_eval_FX];


% 1. STANDARD ADAPTIVE METHOS
tsteps_ST = t_ST(2:end) - t_ST(1:end-1);                    %STEP-SIZE
[~,Y_ST_True] = ode15s(f, t_ST, y0);                        

rel_err_ST = norm(Y_ST_True' - Y_ST) / norm(Y_ST_True');    
Data_ST = [length(t_ST), max(tsteps_ST), min(tsteps_ST), rel_err_ST, num_fun_eval_ST];


% 2.0 MODIFIED ADAPTIVE METHOD VARYING THE FLAG 0:
[~,Y_M0_True] = ode15s(f, t_M0, y0);
    
rel_err_M0 = norm(Y_M0_True' - Y_M0) / norm(Y_M0_True');
Data_M0 = [length(t_M0), max(H0_acp), min(H0_acp), rel_err_M0, num_fun_eval_M0];


% 2.1 MODIFIED ADAPTIVE METHOD VARYING THE FLAG 1:
[~,Y_M1_True] = ode15s(f, t_M1, y0);
    
rel_err_M1 = norm(Y_M1_True' - Y_M1) / norm(Y_M1_True');
Data_M1 = [length(t_M1), max(H1_acp), min(H1_acp), rel_err_M1, num_fun_eval_M1];


% 2.2 MODIFIED ADAPTIVE METHOD VARYING THE FLAG 2:
[~,Y_M2_True] = ode15s(f, t_M2, y0);

rel_err_M2 = norm(Y_M2_True' - Y_M2) / norm(Y_M2_True');
Data_M2 = [length(t_M2), max(H2_acp), min(H2_acp), rel_err_M2, num_fun_eval_M2];


% 2.3 MODIFIED ADAPTIVE METHOD VARYING THE FLAG 3:
[~,Y_M3_True] = ode15s(f, t_M3, y0);
    
rel_err_M3 = norm(Y_M3_True' - Y_M3) / norm(Y_M3_True');
Data_M3 = [length(t_M3), max(H3_acp), min(H3_acp), rel_err_M3, num_fun_eval_M3];


%% PRINT THE RESULTS OF COMPARISON

T = array2table([Data_FX; Data_ST; Data_M0; Data_M1; Data_M2;Data_M3],...
                'VariableNames', ...
                {'Num steps', 'Max step','Min step','Rel_Err','Num fnct Eval'}, ...
                'RowName', ...
                {'Fixed size ERK','Adaptive ERK ST','Adaptive ERK M0','Adaptive ERK M1','Adaptive ERK M2','Adaptive ERK M3'});

fprintf('COMPARE THE METHODS \n')
disp(T)


%% ANALYSIS OF THE MODIFIED ADAPTIVE ERK for flag=0

% 1. PLOTTING SOLUTION
figure(1)
%subplot(3,1,1);
plot(t_M0,Y_M0(1,:),'-*'); hold on;
plot(t_M0,Y_M0(2,:),'-*');

title('SOLUTION ADAPRIVE_ERK MODIFIED');
legend('y1', 'y2')
axis tight;

% 2. PLOTTING THE STEP SIZE
figure(2)
%subplot(3,1,2)
plot(t0_rej,H0_rej,'pentagramr'); hold on;
plot(t_M0(1:end-1),H0_acp,'-ob', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')

title('ACCEPTED/REJECTED SIZE STEP')
legend('Rejected h', 'Accepted h')


% 3. PLOTTING THE ERRORS
exact_loc_err = vecnorm(Y_M0-Y_M0_True') ./ vecnorm(Y_M0_True');    %exact local error

global_err = zeros(length(t_M0),1);                                    %global error
for i=1:length(t_M0)
    global_err(i) = norm(Y_M0(:,1:i) - Y_M0_True(1:i,:)') / norm(Y_M0_True(1:i,:));
end

figure(3)
%subplot(3,1,3)
semilogy(t_M0(2:end), local_err_0, '--o');
hold on;
semilogy(t_M0, exact_loc_err, '--s');
hold on;
semilogy(t_M0, global_err, '--pentagram')

title('ERRORS COMPARISON')
legend('Estimated Local error', 'Exact local error', 'Global error')
