clc
clear

r=1; 
s=1;

oN = [50;300;500];
% oN=50;
omu =0.002:0.002:0.01;
% sigma = 0.01;%1/sqrt(N);     % �Ŷ�ǿ��
oC = [0.1];       % ���Ӷ�
op=0:0.3:0.7;
% op=0.4;

mJ=[];
mJ_DR=[];
mJ_unif=[];
mJ_log=[];
mJ_hn=[];

i=1;
for xh1=1:length(oN)
    N=oN(xh1);
    for xh2=1:length(omu)
        mu=omu(xh2);
        sigma=mu/2;
        for xh3=1:length(oC)
            C=oC(xh3);
            for xh4=1:length(op)
                p=op(xh4);
                [A,A_plus,A_minus,R_plus,R_minus]=ER(N,C,p,mu,sigma,s);
                [ki,k,xm,xi_DR]=DR(R_plus,R_minus,N,s,r,mu,p,C);
                if(k>0.65)
                    continue;
                end
                if(xm<0||sum(xi_DR<=0)>0)
                       continue;
                end
                if(xm>5)
                       continue;
                end
%                 [Rx_mean,xi,F]=x_ode45(A,r,s);
                
                a = min(xi_DR);  % ���ȷֲ�����Сֵ
                b = max(xi_DR);  % ���ȷֲ������ֵ
                xi_unif = a + (b - a) * rand(N, 1);  % ���ɾ��ȷֲ�
                
                mu_log = mean(log(xi_DR(xi_DR > 0)));  % �� xi �������任�����������ֵ
                sigma_log = std(log(xi_DR(xi_DR > 0)));  % �� xi �������任�����������ֵ
                xi_log = lognrnd(mu_log, sigma_log, N, 1);  % ���ɶ�����̬�ֲ�
                
                mu_normal = mean(xi_DR);  % ��̬�ֲ��ľ�ֵ
                sigma_normal = std(xi_DR);  % ��̬�ֲ��ı�׼��
                xi_hn = abs(normrnd(mu_normal, sigma_normal, N, 1));  % ������̬�ֲ��ľ���ֵ
                
%                 J=diag(xi)*A;
                J_DR=diag(xi_DR)*A;
                J_unif=diag(xi_unif)*A;
                J_log=diag(xi_log)*A;
                J_hn=diag(xi_hn)*A;

%                 mJ(i)=max(real(eig(J)));
                mJ_DR(i)=max(real(eig(J_DR)));
                mJ_unif(i)=max(real(eig(J_unif)));
                mJ_log(i)=max(real(eig(J_log)));
                mJ_hn(i)=max(real(eig(J_hn)));
                
                JZ=C*mu*(2*p-1);
                JZ1=C*(mu^2+sigma^2);
                FC=JZ1-JZ^2;
                
                PX=[];
                [PX,~,~]=Tb(xi_DR,C,JZ,FC,s);
                eig_out=sum(xi_DR.*(ki-s))/N;
                yc_mJ_DR(i)=max(max(PX),eig_out);
                
                PX=[];
                [PX,~,~]=Tb(xi_unif,C,JZ,FC,s);
                eig_out=sum(xi_unif.*(ki-s))/N;
                yc_mJ_unif(i)=max(max(PX),eig_out);
                
                PX=[];
                [PX,~,~]=Tb(xi_log,C,JZ,FC,s);
                eig_out=sum(xi_log.*(ki-s))/N;
                yc_mJ_log(i)=max(max(PX),eig_out);
                
                PX=[];
                [PX,~,~]=Tb(xi_hn,C,JZ,FC,s);
                eig_out=sum(xi_hn.*(ki-s))/N;
                yc_mJ_hn(i)=max(max(PX),eig_out);
                
                i=i+1;
            end
        end
    end   
end

i=i-1;

figure
hold on

plot(mJ_DR,yc_mJ_DR,'*')
plot(mJ_unif,yc_mJ_unif,'^')
plot(mJ_log,yc_mJ_log,'o')
plot(mJ_hn,yc_mJ_hn,'s')

x = -1.1:0.1:0.2;
y = x; 
plot(x, y);

xlabel("Solution")
ylabel("SImulaiton")
ans=[mJ_DR;yc_mJ_DR;mJ_unif;yc_mJ_unif;mJ_log;yc_mJ_log;mJ_hn;yc_mJ_hn]'
% save('diff_dist_comp_ER.mat','op','oC','omu','sigma','oN','r','s','mJ_DR','mJ_unif','mJ_log','mJ_hn','yc_mJ_DR','yc_mJ_unif','yc_mJ_log','yc_mJ_hn','i')