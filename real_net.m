clc
clear
s=1;
r=1;
isigma=[0.001,0.003,0.005];
imu=0.005:0.005:0.1;

name=["C0300","C0301","C0306","C0308"];
for i=1:4
    m=name(i);
    load(m+'.mat')
    [N,~]=size(adj_matrix);
    A_plus=(adj_matrix>0);
    A_plus=A_plus+A_plus';
    A_minus=(adj_matrix<0);
    A_minus=A_minus+A_minus';
    C=sum(sum(A_plus+A_minus))/N/(N-1);
    p=sum(sum(A_plus))/sum(sum(A_plus+A_minus));
    S_J=[];
    S_PX=[];
    for xh1=1:length(isigma)
        sigma=isigma(xh1);
        R=sigma.*randn(N,N);
        for xh2=1:length(imu)
            mu=imu(xh2);
            R_plus=A_plus.*(R+mu);
            R_minus=A_minus.*(R+mu);
            A=R_plus-R_minus-s*eye(N);
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
    save('real_net_'+m+'.mat','imu','isigma','S_J','S_PX')
    figure
    hold on
    for xh1=1:length(isigma)
        len=length(S_J(xh1,:));
        plot(imu(1:len),S_J(xh1,:),'o')
        plot(imu(1:len),S_PX(xh1,:),'-')
 
    end
    hold off
end

