%%
clc
clear
s=1;
r=1;
C=0.1;
N=50;
mu=0.01;
sigma=0.01;
p=1;

max_ki = Inf;
while max_ki >= s
    [A,A_plus,A_minus,R_plus,R_minus]=ER(N,C,p,mu,sigma,s);
    [ki,k,xm,xi]=DR(R_plus,R_minus,N,s,r,mu,p,C);
    max_ki=max(abs(ki));
end
max(ki)
min(ki)
[Rx_mean,Rxi,F]=x_ode45(A,r,s);

[i,j]=size(F);
figure;
hold on
for ii=1:j
   plot(F(:,ii),'color',[189 190 192]/255,'linewidth',1); 
   plot(i,xi(ii),'r.','color',[239 124 120]/255,'linewidth',1,'MarkerSize',8);
end

J=diag(xi)*A;
RJ=diag(Rxi)*A;
ej=eig(J);
rJ=real(ej);
iJ=imag(ej);

JZ=C*mu*(2*p-1);
JZ1=C*(mu^2+sigma^2);
FC=JZ1-JZ^2;
[PX,PY,T_values]=Tb(xi,C, JZ, FC,s);

eig_out=sum(xi.*(ki-s))/N;

figure;
hold on
% plot(real(eig(RJ)),imag(eig(RJ)),'r*')
plot(rJ,iJ,'b*')
plot(PX,PY,'ro')
plot(-xi*(s+JZ),0,'yo')
plot(eig_out,0,'^')
F_xi=-xi;
% save('eig_out.mat','C','iJ','rJ','mu','sigma','p','N','xi','F_xi','PX','PY','eig_out')
% save('eig_xi.mat','C','iJ','rJ','mu','sigma','p','N','xi','F_xi','PX','PY','eig_out')