# Two Dimensional Wavelet Based Statistical Monitoring of Image Data
#### Mehdi Koosha, Rassoul Noorossana, Orod Ahmadi

### Abstract
Advances in imaging technology has paved the grounds to capture and analyze product images from different perspectives in real time. In high volume and highly automated processes, effective methods are required to monitor several quality characteristics of interest. Statistical process monitoring techniques are widely used methods to monitor the stability of processes. Image data could be used for statistical process monitoring. The advantages of applying images in monitoring processes, have led to their increasing use. When images are used, complicated structures of data need to be analyzed. Various methods including data transformation and dimension reduction may be used to cope with these complexities. In this paper, a two dimensional wavelet is applied to images. Using this transformation estimated wavelet coefficients are obtained. To monitor the expected values of these coefficient, a likelihood ratio test is developed. The proposed approach is able to detect out-of-control conditions and estimate change points effectively. This diagnostic information can help practitioners to find the assignable cause(s) of the change and return the process back to the in-control state. One of the advantages of the proposed method compared to other competing methods is that the performance of the method does not depend on the fault location. Simulation studies indicate satisfactory performance of proposed method for detecting out-of-control conditions and change point estimation. The performance of the proposed method is also compared with two existing approach and results indicate superior performance for the proposed method.

## Deployment

### 1- Preprocessing code

#### a- Reading initial image
```
  I=imread('Tile_Image_00000.jpg');% reading the initial image
```
The initial image can be accessed in [GitHub](https://raw.githubusercontent.com/mehdikoosha/SPC-image-data-wavelets/master/Tile_Image_00000.JPG)


#### b- Change image type to grayscale
```
I=rgb2gray(I);%change the initial image to grayscale
```
#### c- Image preprocessing
```
pixs=[256;256]; %size of the output image
I=I(160:2300, 740:2880);% Omitting the useless parts of the image borders
I=imresize(I,pixs');%resize the image to the preferred size
BG=imopen(I,strel('disk',15)); %Estimating the Value of Background Pixels
Iunfrm=imsubtract(I,BG); %Create an Image with a Uniform Background
Iadjust=imadjust(Iunfrm); %Adjusting the contrast-see imadjust for details
Ires=imresize(Iadjust,pixs'); %Resizing
I=Ires;

```
The resulted image is shown in here:
![](https://raw.githubusercontent.com/mehdikoosha/SPC-image-data-wavelets/master/Tile_Nom.jpg)

### 2- 2Dwavelet_Phase 1 code
#### a- Setting initial parameters
```
n=256; %size of image   
SNR=25; %determined signal to noise ratio
n_level=2; %level of decomposition
UCL=4.1; %Upper control Limit
name=haar; %wavelet function type
```
#### b- Estiamting image standard deviation
```
I=imread('Tile_Nom.jpg');
I=im2double(I);
Signal=std(I(:));
sigma_error=Signal/(10^(SNR/20));
```
Sigma error is computed eqaul to 0.0118 in our case.
#### c- Phase I parameter estimation
In this step, Phase I images are generated and Mean vector and standard deviation vector are estimated
```
m=10000; %number of Phase I in-control images
for a=1:m
    I0=imnoise(I,'gaussian',0,sigma_error^2);
    tetha=wavedec2(I0,n_level,name);
    tetha(a,:)=tetha(1:n^2);
end;
tetha0=mean(tetha);
sigma_0=std(tetha);
```
#### d- ARL estimation based on the predefined UCL
In this step, the in control average run length (ARL) is evaluated based on the predefined UCL in step 2a. The UCL is changed until the predetermined in control ARL (in this paper, 200) is achieved.
```
for m=1:1000
    clear LAMDA Lamda
    Loop=0;
    T=0;
    while Loop==0
        T=T+1
        I01=imnoise(I,'gaussian',0,sigma_error^2);
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
                Lamda(T,tau)=(((tau-T)/2)*log(1+2*GAMMA2(T,tau)/(n^2)))+...
                ...(2*GAMMA2(T,tau)*sum(w_p(tau+1:T).^2)+...
                ...n*GAMMA2(T,tau)*sqrt(2)*sum(w_p(tau+1:T))-...
                ...GAMMA2(T,tau)^2*(T-tau)/2)/(2*n^2+4*GAMMA2(T,tau));
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
```
The Upper Control Limit to achieve in-control ARL of approximately 200 is computed as 4.1 based on trial and error. The in-control ARL is 200.98 in this circumstance.
### 3- 2Dwavelet_phase II code
In this step, shifts in multiple sizes, magnitudes and in diffrent locations are implemented in the data and the results are reported for each setting.
```
n=256; %size of image   
SNR=25; %determined signal to noise ratio
n_level=2; %level of decomposition
UCL=4.1; %Upper control Limit
name=haar; %wavelet function type
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
MEDDEV=median(DEV);
RESULTS(Z,:)=[ARL,STDRL,MEDRL,ADEV,STDDEV,MEDDEV];
end
```
Results can be found in Tables of Section 4 of the paper. It should be noted that comparison is done with two competing methods including Koosha et al. (2017) and Amirkhani and Amiri (2019). The codes of the competing methods can be accessed at the following links:

Koosha et al. (2017): 
https://github.com/mehdikoosha/SPC-image-data-wavelets

Amirkhani and Amiri (2019): https://github.com/fzdamirkhani/QREI-2019

The permanenet DOI for this repository is [![DOI](https://zenodo.org/badge/355555897.svg)](https://zenodo.org/badge/latestdoi/355555897)

