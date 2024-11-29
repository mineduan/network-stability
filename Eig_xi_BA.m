clc
clear

r=1; 
s=1;

oN = [50;100;300;500];
omu =0.002:0.002:0.02;
%1/sqrt(N);     % 扰动强度
% oC = [0.1,0.3,0.5,0.7,1];       % 连接度
op=0:0.2:1;
S=[];
RS=[];
i=1;
for xh1=1:length(oN)
    N=oN(xh1);
    m0=round(N/2);
    m=round(N/4);
    for xh2=1:length(omu)
        mu=omu(xh2);
        sigma = mu/2;
%         for xh3=1:length(oC)
%             C=oC(xh3);
            for xh4=1:length(op)
                p=op(xh4);
                [A,A_plus,A_minus,R_plus,R_minus,C,p]=BA(m0,m,N,p,mu,sigma);
                [ki,k,xm,xi]=DR(R_plus,R_minus,N,s,r,mu,p,C);
                if(abs(k)>0.65)
                    continue;
                end
                if(max(abs(ki))>1)
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
%                 S(i)=-(s+JZ)*min(xi);
                
                i=i+1;
                
                
            end
%         end
    end   
end

figure
hold on
plot(S,RS,'.');

x = -1.1:0.1:0.2;
y = x; 
plot(x, y);

xlabel("Solution")
ylabel("SImulaiton")
% save('Eig_xi_BA_comp.mat','op','omu','sigma','oN','r','s','S','RS','i')