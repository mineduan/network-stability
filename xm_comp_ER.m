clc
clear
close all

s=1;
C=0.1;

r=1;
N=300;
mu=0.04;
sigma=0.01;
p=0:0.01:1;
%p=[p,p,p];
for i=1:length(p)
    [A,A_plus,A_minus,R_plus,R_minus]=ER(N,C,p(i),mu1,mu2,sigma1,sigma2,s);
    [xm,xi]=DR(R_plus,R_minus,N,s,r,mu1,mu2,p(i),C);
    [Rxm,~,~]=x_ode45(A,r,s);

    Axm(i)=xm;
    ARxm(i)=Rxm;
end

% figure for xm_comp on p
figure(1)
hold on
plot(p,Axm,'o');
plot(p,ARxm,'-');
legend('simulation','solution');
set(gca,'FontSize',16);
set(gca,'FontName','Times New Roman');
xlabel('p','FontSize',16, 'FontWeight', 'bold','Fontangle','italic')
ylabel('<x>','FontSize',16, 'FontWeight', 'bold','Fontangle','italic')

hold off
% save('xm_comp_for_p.mat','C','N','mu1','p','sigma1','T1')
%%
clc
clear
close all

s=1;
C=0.1;

r=1;
mu1=0.01;
mu2=0.01;

N=100:10:1000;
p=0.3;
xh=1;
%p=[p,p,p];
for i=1:length(N)
    sigma1=1/sqrt(N(i));
    sigma2=1/sqrt(N(i));
    Sxm=0;
    SRxm=0;
    for j=1:xh
            [A,A_plus,A_minus,R_plus,R_minus]=ER(N(i),C,p,mu1,mu2,sigma1,sigma2,s);
            [xm,xi]=DR(R_plus,R_minus,N(i),s,r,mu1,mu2,p,C);
            [Rxm,~,~]=x_ode45(A,r,s);
            Sxm=Sxm+xm;
            SRxm=SRxm+Rxm;
    end
    Axm(i)=Sxm/xh;
    ARxm(i)=SRxm/xh;
end


T1=[N;Axm;ARxm];
% T2=sortrows(T1',1);
% T1=T2';

%% figure
figure(1)
hold on
plot(T1(1,:),T1(3,:),'^');
plot(T1(1,:),T1(2,:),'-');
legend('simulation','solution');
set(gca,'FontSize',16);
set(gca,'FontName','Times New Roman');
xlabel('p','FontSize',16, 'FontWeight', 'bold','Fontangle','italic')
ylabel('<x>','FontSize',16, 'FontWeight', 'bold','Fontangle','italic')
hold off
% save('xm_comp_for_N.mat','C','N','mu1','p','sigma1','T1')