clear
clc
%% parameters
n=256; %size of image   
SNR=25; %determined signal to noise ratio
n_level=2; %level of decomposition
UCL=4.1; %Upper control Limit
name=haar; %wavelet function type
%% reading initial iamge
I=imread('Tile_Nom.jpg');
I=im2double(I);
Signal=std(I(:));
sigma_error=Signal/(10^(SNR/20));
%% Phase I
m=10000;
for a=1:m
    I0=imnoise(I,'gaussian',0,sigma_error^2);
    tetha=wavedec2(I0,n_level,name);
    tetha(a,:)=tetha(1:4096);
end;
tetha0=mean(tetha);
sigma_0=std(tetha);
for m=1:1000
clear LAMDA Lamda
Loop=0;
T=0;
while Loop==0
    T=T+1
    I01=imnoise(I,'gaussian',0,sigma_error^2);
    tethA_2(T,:)=wavedec2(I01,n_level,name);
    TETHA_2=tethA_2(T,:);
    tetha_2(T,:)= TETHA_2(1:4096);  
    w(T)=sum(((tetha_2(T,:)-tetha0)./sigma_0).^2);
    w_p(T)=(w(T)-n^2)/(sqrt(2)*n);
    if T>1
        for tau=1:T-1
            tetha1=mean((tetha_2(tau+1:T,:)),1);
            GAMMA2(T,tau)=sum(((tetha1-tetha0)./sigma_0).^2);
            GAMMA2(T,tau)=GAMMA2(T,tau)-(n^2/(T-tau));
            Lamda(T,tau)=(((tau-T)/2)*log(1+2*GAMMA2(T,tau)/(n^2)))+(2*GAMMA2(T,tau)*sum(w_p(tau+1:T).^2)+n*GAMMA2(T,tau)*sqrt(2)*sum(w_p(tau+1:T))-GAMMA2(T,tau)^2*(T-tau)/2)/(2*n^2+4*GAMMA2(T,tau));
        end
        [LAMDA(T),index]=max(Lamda(T,:));
        if LAMDA(T)>UCL
            Loop=1;
            RL=T-Change_Point;
            RL(m)=RL;
        end
    end
end
end
ARL=mean(RL)
