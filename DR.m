function [ki,k,xm,xi]=DR(R_plus,R_minus,N,s,r,mu,p,C)
    
%     k_plus=sum(R_plus,1); %求度行和
%     k_minus=sum(R_minus,1);
    ki=sum(R_plus,2)-sum(R_minus,2);
    k=mean(ki);

%     xm=r/(s+(N-1)*C*mu*(1-2*p));
    xm=r/(s-k);
    xi=(r+ki*xm)/s;
    xi(xi<0)=0;
end