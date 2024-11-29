clc
clear

name=["C0300","C0301","C0306","C0308"];
for i=1:4
    m=name(i);
   load('real_net_'+m+'.mat') 
   figure
    hold on
    for xh1=1:length(isigma)
        len=length(S_J(xh1,:));
        plot(imu(1:len),S_J(xh1,:),'o')
        plot(imu(1:len),S_PX(xh1,:),'-','LineWidth',2)
    end
    hold off
    if i==1
        F1=[];
        for xh1=1:length(isigma)
            len=length(S_J(xh1,:));
            F1=[F1,imu(1:len)',S_J(xh1,:)',S_PX(xh1,:)'];
        end
    elseif i==2
        F2=[];
        for xh1=1:length(isigma)
            len=length(S_J(xh1,:));
            F2=[F2,imu(1:len)',S_J(xh1,:)',S_PX(xh1,:)'];
        end
        
    elseif i==3
        F3=[];
        for xh1=1:length(isigma)
            len=length(S_J(xh1,:));
            F3=[F3,imu(1:len)',S_J(xh1,:)',S_PX(xh1,:)'];
        end
        
    elseif i==4
        F4=[];
        for xh1=1:length(isigma)
            len=length(S_J(xh1,:));
            F4=[F4,imu(1:len)',S_J(xh1,:)',S_PX(xh1,:)'];
        end
    end
        
end