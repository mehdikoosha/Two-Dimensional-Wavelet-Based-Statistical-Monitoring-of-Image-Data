clear;
clc;
%% parameters
n=64;
SNR=25;
n_level=2;
Change_Point=5;
UCL=8.1;
name='haar';
%% reading initial iamge
I=imread('Tile_Nom.jpg');
I=im2double(I);
Signal=std(I(:));
%sigma_error=0.01;
sigma_error=Signal/(10^(SNR/20));
%%Phase I
m=10000;
for a=1:m/10
    I0=imnoise(I,'gaussian',0,sigma_error^2);
    tetha=wavedec2(I0,n_level,name);
    tethaA(a,:)=tetha(1:4096);
    a
end;
for b=1:m/10
    I0=imnoise(I,'gaussian',0,sigma_error^2);
    tetha=wavedec2(I0,n_level,name);
    tethaB(b,:)=tetha(1:4096);
    b
end;
for c=1:m/10
    I0=imnoise(I,'gaussian',0,sigma_error^2);
    tetha=wavedec2(I0,n_level,name);
    tethaC(c,:)=tetha(1:4096);
    c
end;
for d=1:m/10
    I0=imnoise(I,'gaussian',0,sigma_error^2);
    tetha=wavedec2(I0,n_level,name);
    tethaD(d,:)=tetha(1:4096);
    d
end;
for e=1:m/10
    I0=imnoise(I,'gaussian',0,sigma_error^2);
    tetha=wavedec2(I0,n_level,name);
    tethaE(e,:)=tetha(1:4096);
    e
end;
for f=1:m/10
    I0=imnoise(I,'gaussian',0,sigma_error^2);
    tetha=wavedec2(I0,n_level,name);
    tethaF(f,:)=tetha(1:4096);
    f
end;
for g=1:m/10
    I0=imnoise(I,'gaussian',0,sigma_error^2);
    tetha=wavedec2(I0,n_level,name);
    tethaG(g,:)=tetha(1:4096);
    g
end;
for h=1:m/10
    I0=imnoise(I,'gaussian',0,sigma_error^2);
    tetha=wavedec2(I0,n_level,name);
    tethaH(h,:)=tetha(1:4096);
    h
end;
for i=1:m/10
    I0=imnoise(I,'gaussian',0,sigma_error^2);
    tetha=wavedec2(I0,n_level,name);
    tethaI(i,:)=tetha(1:4096);
    i
end;
for j=1:m/10
    I0=imnoise(I,'gaussian',0,sigma_error^2);
    tetha=wavedec2(I0,n_level,name);
    tethaJ(j,:)=tetha(1:4096);
    j
end;
tetha=[tethaA;tethaB;tethaC;tethaD;tethaE;tethaF;tethaG;tethaH;tethaI;tethaJ];
tetha0=mean(tetha);
sigma_0=std(tetha);
%% Phase II
for m=1:1000
clear LAMDA Lamda
Loop=0;
T=0;
while Loop==0
    T=T+1
    I01=imnoise(I,'gaussian',0,sigma_error^2);
    if T>Change_Point
        Shift=zeros(256,256);
        delta=4/255;
        % shift magnitude
        %%%%% shift location
        a=10;
        b=10;
        llength=10;
        for i=a:a+llength
            for j=b:b+llength
                Shift(i,j)=delta;
            end
        end
        I01=I01+Shift;
    end
    tethA_2(T,:)=wavedec2(I01,n_level,name);
    TETHA_2=tethA_2(T,:);
    tetha_2(T,:)= TETHA_2(1:4096);  
    w(T)=sum(((tetha_2(T,:)-tetha0)./sigma_0).^2);
    w_p(T)=(w(T)-n^2)/(sqrt(2)*n);
    if T>1
        for tou=1:T-1
            tetha1=mean((tetha_2(tou+1:T,:)),1);
            GAMMA2(T,tou)=sum(((tetha1-tetha0)./sigma_0).^2);
            GAMMA2(T,tou)=GAMMA2(T,tou)-(n^2/(T-tou));
            Lamda(T,tou)=(((tou-T)/2)*log(1+2*GAMMA2(T,tou)/(n^2)))+(2*GAMMA2(T,tou)*sum(w_p(tou+1:T).^2)+n*GAMMA2(T,tou)*sqrt(2)*sum(w_p(tou+1:T))-GAMMA2(T,tou)^2*(T-tou)/2)/(2*n^2+2*GAMMA2(T,tou));
        end
        [LAMDA(T),index]=max(Lamda(T,:));
        if LAMDA(T)>UCL
            Loop=1;
            index;
            DEV=index-Change_Point-1;
            RL=T-Change_Point-1;
            RLL(m)=RL;
            DEVV(m)=max(0,DEV);
        end
    end
end
end
ARL=mean(RLL)
STDRL=std(RLL)
ADEV=mean(DEVV)