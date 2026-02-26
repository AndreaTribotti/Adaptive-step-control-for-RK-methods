function [t,y,h] = IRK(f,tspan,y0,n,butcher)
b=butcher.b; A=butcher.A; c=butcher.c;
q=length(y0); % serve per definire bene F e il for orribile
l=length(b); % serve per definire bene il for orribile
ql=q*l; % serve in NWT_solver_reshape (ho finito le lettere)
y(:,1)=y0; t(1)=tspan(1);
h=(tspan(2)-tspan(1))/n;
for k=2:n
    tt=t(k-1);
    F = @(x) cell2mat(arrayfun(@(i) ...
    x((i-1)*q+1:i*q) - f(tt + h * c(i), ...
    y(:,k-1) + h * sum(A(i, :) .* reshape(x, q, l), 2)), ...
    1:l, 'UniformOutput', false));
    [K,N]=NWT_solver_reshape(rand(ql,1),F,10^-3,1000);
    yk=y(:,k-1);
    % size(yk)
    % size(K(1+3*(1-1):3*1))
    for s=1:l
        yk=yk+h*b(s)*K(1+q*(s-1):s*q); % questo for è orribile, sarebbe da cambiare
    end
    y(:,k)=yk;
    t(k)=t(k-1)+h;

end