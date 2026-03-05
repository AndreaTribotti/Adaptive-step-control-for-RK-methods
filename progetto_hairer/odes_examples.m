syms t y(t)

% f1
ode1= diff(y,t)==3+cos(t)-y; cond1= y(0)==1;
ysolv1(t)=dsolve(ode1,cond1);
sol1=matlabFunction(ysolv1);
f1=@(t,y)(3+cos(t)-y); 
% figure
% fplot(ysolv1,[0 20])
% title("f1")


% f2
ode2= diff(y,t)==y*(2+y)*(2-y); cond2= y(0)==5;
ysolv2(t)=dsolve(ode2,cond2);
f2=@(t,y)(y.*(2+y).*(2-y)); 
% figure
% fplot(ysolv2,[0 2])
% title("f2")