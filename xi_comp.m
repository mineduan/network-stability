%%  ER
clc
clear

N=50;
mu=0.1;
sigma=0.01;
p=0;
s=1;
r=1;
C=0.1;
[A,A_plus,A_minus,R_plus,R_minus]=ER(N,C,p,mu,sigma,s);
[ki,k,xm,xi]=DR(R_plus,R_minus,N,s,r,mu,p,C);
max(abs(ki))
[Rx_mean,Rxi,F]=x_ode45(A,r,s);
sum(Rxi<10e-4)
sum(xi<=0)
[i,j]=size(F);
figure;
hold on
for ii=1:j
   plot(F(:,ii),'color',[189 190 192]/255,'linewidth',1); 
   plot(i,xi(ii),'r.','color',[239 124 120]/255,'linewidth',1,'MarkerSize',8);
end
set(gca,'FontSize',16);
set(gca,'FontName','Times New Roman');
xlabel('T','FontSize',16, 'FontWeight', 'bold','Fontangle','italic')
ylabel('x_i','FontSize',16, 'FontWeight', 'bold','Fontangle','italic')
% legend('simulation','solution')
hold off

%% BA
clc
clear

N=50;
mu=0.03;
sigma=0;

p=0;
s=1;
r=1;
m0=25;
m=5;
[A,A_plus,A_minus,R_plus,R_minus,C,p]=BA(m0,m,N,p,mu,sigma);
[ki,k,xm,xi]=DR(R_plus,R_minus,N,s,r,mu,p,C);
max(ki)
min(ki)
[Rx_mean,Rxi,F]=x_ode45(A,r,s);

sum(Rxi<10e-4)
sum(xi<=0)
[i,j]=size(F);

% figure
figure;
hold on
for ii=1:j
   plot(F(:,ii),'color',[189 190 192]/255,'linewidth',1); 
   plot(i,xi(ii),'r.','color',[239 124 120]/255,'linewidth',1,'MarkerSize',8);
end
set(gca,'FontSize',16);
set(gca,'FontName','Times New Roman');
xlabel('T','FontSize',16, 'FontWeight', 'bold','Fontangle','italic')
ylabel('x_i','FontSize',16, 'FontWeight', 'bold','Fontangle','italic')
% legend('simulation','solution')
hold off
% save('xi_comp_BA_N50.mat','N','p','m0','m','C','A','xm','xi','F')