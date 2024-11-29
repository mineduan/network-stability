clc
clear

r=1; 
s=1;

oN = [50;100;300;500];
omu =0.002:0.002:0.02;
sigma = 0.01;%1/sqrt(N);     % 扰动强度
oC = [0.1,0.3,0.5];       % 连接度
op=0:0.2:1;
S=[];
RS=[];
i=1;
for xh1=1:length(oN)
    N=oN(xh1);
    for xh2=1:length(omu)
        mu=omu(xh2);
        for xh3=1:length(oC)
            C=oC(xh3);
            for xh4=1:length(op)
                p=op(xh4);
                [A,A_plus,A_minus,R_plus,R_minus]=ER(N,C,p,mu,sigma,s);
                [ki,k,xm,xi]=DR(R_plus,R_minus,N,s,r,mu,p,C);
                if(k>0.65)
                    continue;
                end
                if(xm<0||sum(xi<=0)>0)
                       continue;
                end
                if(xm>5)
                       continue;
                end
                J=diag(xi)*A;
                RS(i)=max(real(eig(J)));
                
                JZ=C*mu*(2*p-1);
                
                eig_out=sum(xi.*(ki-s))/N;
                
                S(i)=max(-(s+JZ)*min(xi),eig_out);

                i=i+1;
            end
        end
    end   
end

figure
hold on
plot(S,RS,'.');

x = -1.1:0.1:0.2;
y = x; 
plot(x, y);

xlabel("Solution")
ylabel("Simulaiton")
% save('Eig_xi_ER_comp.mat','op','oC','omu','sigma','oN','r','s','S','RS','i')