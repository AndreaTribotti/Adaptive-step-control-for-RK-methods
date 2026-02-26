function [t,u,num_fval,b,c,A] = ERK3_2inputs(f,tspan,y0,h,b1,b2)
u(:,1)=y0; t(1)=tspan(1);
N=round((tspan(2)-tspan(1))/h); 

% Trovo tutti i coefficienti dati b1 e b2
b3 = 1 - b1 - b2; 
fun = @(x) [
    b2 * x(1) + b3 * x(2) - 1/2;           % Eq 2: b2*c2 + b3*c3 = 1/2
    b2 * x(1)^2 + b3 * x(2)^2 - 1/3;       % Eq 3: b2*c2^2 + b3*c3^2 = 1/3
    x(2) * x(5) * b3 - 1/6;                % Eq 4: c3 * a32 * b3 = 1/6
    x(1) - x(3);                           % Eq 5: c2 = a21
    x(2) - ( x(4) + x(5))                  % Eq 6: c3 = a31 + a32
];
options = optimoptions('fsolve', 'Display', 'none');
x0 = [0.1; 0.1; 0.1; 0.1; 0.1];
[x, ~, ~] = fsolve(fun, x0, options);
c2 = x(1);
c3 = x(2);
a21 = x(3);
a31 = x(4);
a32 = x(5);

% Creo b,c,A
b=[b1,b2,b3];
c=[0,c2,c3];
A(2,1)=a21; A(3,1)=a31; A(3,2)=a32;

for k=2:N+1
    K1=f(      t(k-1),               u(:,k-1) );
    K2=f(   t(k-1)+c(2)*h,       u(:,k-1)+h*A(2,1)*K1 );
    K3=f(   t(k-1)+c(3)*h,   u(:,k-1)+h*A(3,1)*K1+h*A(3,2)*K2 );
    u(:,k)=u(:,k-1)+h*(b(1)*K1+b(2)*K2+b(3)*K3);
    t(k,1)=t(k-1)+h;
end
num_fval=3*N;
end