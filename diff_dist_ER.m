clc
clear

r=1; 
s=1;

oN = [50;300;500];
omu =0.002:0.002:0.01;
% sigma = 0.01;%1/sqrt(N);     % 扰动强度
oC = [0.1];       % 连接度
op=0:0.4:1;
S=[];
RS=[];
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
                [Rx_mean,xi,F]=x_ode45(A,r,s);
                
                a = min(xi_DR);  % 均匀分布的最小值
                b = max(xi_DR);  % 均匀分布的最大值
                xi_unif = a + (b - a) * rand(1, N);  % 生成均匀分布
                
                mu_log = mean(log(xi_DR(xi_DR > 0)));  % 对 xi 做对数变换，避免对数负值
                sigma_log = std(log(xi_DR(xi_DR > 0)));  % 对 xi 做对数变换，避免对数负值
                xi_log = lognrnd(mu_log, sigma_log, 1, N);  % 生成对数正态分布
                
                mu_normal = mean(xi_DR);  % 正态分布的均值
                sigma_normal = std(xi_DR);  % 正态分布的标准差
                xi_hn = abs(normrnd(mu_normal, sigma_normal, 1, N));  % 生成正态分布的绝对值
                
                J=diag(xi)*A;
                J_DR=diag(xi_DR)*A;
                J_unif=diag(xi_unif)*A;
                J_log=diag(xi_log)*A;
                J_hn=diag(xi_hn)*A;

                mJ(i)=max(real(eig(J)));
                mJ_DR(i)=max(real(eig(J_DR)));
                mJ_unif(i)=max(real(eig(J_unif)));
                mJ_log(i)=max(real(eig(J_log)));
                mJ_hn(i)=max(real(eig(J_hn)));
                
                i=i+1;
            end
        end
    end   
end

i=i-1;

figure
hold on

plot(mJ,mJ_DR,'*')
plot(mJ,mJ_unif,'^')
plot(mJ,mJ_log,'o')
plot(mJ,mJ_hn,'s')

sum(abs(mJ-mJ_DR))/i
sum(abs(mJ-mJ_unif))/i
sum(abs(mJ-mJ_log))/i
sum(abs(mJ-mJ_hn))/i

x = -1.1:0.1:0.2;
y = x; 
plot(x, y);

xlabel("Solution")
ylabel("SImulaiton")
res=[mJ;mJ_DR;mJ_unif;mJ_log;mJ_hn]'
% save('diff_dist_ER.mat','op','oC','omu','sigma','oN','r','s','mJ','mJ_DR','mJ_unif','mJ_log','mJ_hn','i')