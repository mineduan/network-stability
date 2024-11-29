clc
clear
s=1;
r=1;
% iN=[100,300,500,800];
N=300;
% ip=0:0.05:1;
p=0;
iC=[0.01,0.1,0.3];
% mu=0.04;
imu=0.002:0.002:0.05;
sigma=0.005;


S_J=[];
S_PX=[];
for xh1=1:length(iC)
    C=iC(xh1);
    for xh2=1:length(imu)
        mu=imu(xh2);
        [A,A_plus,A_minus,R_plus,R_minus]=ER(N,C,p,mu,sigma,s);
        [ki,k,xm,xi]=DR(R_plus,R_minus,N,s,r,mu,p,C);
        
        if(xm<0||sum(xi<=0)>0)
            break;
        end
        if(xm>5)
            break;
        end
        
        J=diag(xi)*A;
        S_J(xh1,xh2)=max(real(eig(J)));
        
        JZ=C*mu*(2*p-1);
        JZ1=C*(mu^2+sigma^2);
        FC=JZ1-JZ^2;
        [PX,PY,T_values]=Tb(xi, C, JZ, FC,s);
        eig_out=sum(xi.*(ki-s))/N;
        S_PX(xh1,xh2)=max(max(PX),eig_out);
        
        %             figure;
        %             hold on
        %             plot(real(eig(J)),imag(eig(J)),'b*')
        %             plot(PX,PY,'ro')
        %             plot(eig_out,0,'^')
    end
end

figure
hold on
for xh1=1:length(iC)
    len=length(S_J(xh1,:));
    plot(imu(1:len),S_J(xh1,:),'o')
    plot(imu(1:len),S_PX(xh1,:),'-')
    
end
legend('N=100','N=300','N=500','N=800')
hold off
% save('S_C_mu_1.mat','iC','imu','sigma','N','p','S_J','S_PX')

F=[];
for xh1=1:length(iC)
    len=length(S_J(xh1,:));
    F=[F,imu(1:len)',S_J(xh1,:)',S_PX(xh1,:)'];
end