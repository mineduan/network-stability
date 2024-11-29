function [x_mean,xi,F] = x_ode45(A,r,s)
x_initial=0.5*ones(1,length(A(1,:)));
options=odeset('Reltol',1e-5);
[ij,F]=ode45(@funX,[0,100],x_initial,options,A,r,s);
xi=F(end,:);
x_mean=sum(xi)/length(xi);