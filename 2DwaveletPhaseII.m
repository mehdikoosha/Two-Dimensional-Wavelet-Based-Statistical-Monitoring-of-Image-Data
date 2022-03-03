clear;
clc;
%% parameters
n=256;
SNR=25;
n_level=2;
Change_Point=5;
UCL=4.1;
name='haar';
%% reading initial iamge
I=imread('Tile_Nom.jpg');
I=im2double(I);
Signal=std(I(:));
%sigma_error=0.01;
sigma_error=Signal/(10^(SNR/20));
%%Phase I
m=10000;
for a=1:m
    I0=imnoise(I,'gaussian',0,sigma_error^2);
    tetha=wavedec2(I0,n_level,name);
    tethaA(a,:)=tetha(1:n^2);
    a
end;
tetha0=mean(tetha);
sigma_0=std(tetha);
%% Phase II
for m=1:1000
clear LAMDA Lamda
Loop=0;
T=0;
for Z=[1:10]
while Loop==0
    T=T+1
    I01=imnoise(I,'gaussian',0,sigma_error^2);
    if T>Change_Point
        Shift=zeros(pixs(1),pixs(2));
        delta=Z/255; % shift magnitude
        %%%%% shift location
        Shift_horizontal_location=[10,20];
        Shift_vertical_Location=[10,20];
        for i= Shift_horizontal_location(1): Shift_horizontal_location(2)
            for j=Shift_vertical_Location(1):Shift_vertical_Location(2)
                Shift(i,j)=delta;
            end
        end
        I01=I01+Shift;
        for i=1:pixs(1)
            for j=1:pixs(2)
                I01(i,j)=max(I01(i,j),0); %To guarantee that each pixel's intensity is not smaller than zero
                I01(i,j)=min(I01(i,j),255); %To guarantee that each pixel's intensity is not larger than 255
            end
        end
    end    
    tethA_2(T,:)=wavedec2(I01,n_level,name);
    TETHA_2=tethA_2(T,:);
    tetha_2(T,:)= TETHA_2(1:n^2);  
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
            index;
            RL=T-Change_Point;
            DEV=abs(index-Change_Point);
            RL(m)=RL;
            DEV(m)=DEV;
        end
    end
end
ARL=mean(RL);
STDRL=std(RL);
MEDRL=median(RL);
ADEV=mean(DEV);
STDDEV=std(DEV);
MEDDEV=median(dev);
RESULTS(Z,:)=[ARL,STDRL,MEDRL,ADEV,STDDEV,MEDDEV];
end