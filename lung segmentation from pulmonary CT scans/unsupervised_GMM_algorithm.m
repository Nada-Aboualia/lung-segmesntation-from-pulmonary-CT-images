a1=double(imread('CT_001.bmp'));
a2=double(imread('CT_002.bmp'));
a3=double(imread('CT_003.bmp'));
a4=double(imread('CT_004.bmp'));

g1=double(imread('Ground_Truth_CT_001.bmp'));
g2=double(imread('Ground_Truth_CT_002.bmp'));
g3=double(imread('Ground_Truth_CT_003.bmp'));
g4=double(imread('Ground_Truth_CT_004.bmp'));

% Initialization

init_mu1 = 30;
init_segma = 100;
init_prob = 0.5;
init_mu2 = 200;
segma2_old = 100;
Prob2_old = 0.5;

A = [a1 a2 a3 a4];
x= [g1 g2 g3 g4];
[h ,w] = size(A);


% Reading and Calculating the Sizefor index = 1:50
for index=1:50       % Number of Iterations
for k1 = 1:h
for k2 = 1:w
gray = A(k1,k2);
if gray < 255
P1 = init_prob.*Gauss_aym(gray,init_mu1,init_segma);
P2 = Prob2_old.*Gauss_aym(gray,init_mu2,segma2_old);
pi_1(k1,k2) = P1./(P1 + P2);
pi_2(k1,k2) = P2./(P1 + P2);
else
pi_1(k1,k2) =0;
pi_2(k1,k2) = 0;
end
end
end
  
% Calculate the responsibilityProb_new 
Prob1_new = sum(sum(pi_1))./(sum(sum(pi_1)) + sum(sum(pi_2)));
Prob2_new = sum(sum(pi_2))./(sum(sum(pi_1)) + sum(sum(pi_2)));
%claculate new mean and variance
%Mean new values
Mu1_new = (sum(sum(pi_1.*A)))./(sum(sum(pi_1)));
Mu2_new = (sum(sum(pi_2.*A)))./(sum(sum(pi_2)));
%Variance new values
Var1_new = (sum(sum(pi_1.*((A-init_mu1).^2))))./(sum(sum(pi_1)));
Var2_new = (sum(sum(pi_2.*((A-init_mu2).^2))))./(sum(sum(pi_2)));
init_mu1 = Mu1_new;
init_segma = Var1_new;
init_prob = Prob1_new;
init_mu2 = Mu2_new;
segma2_old = Var2_new;
Prob2_old = Prob2_new;
Mu_1_Vec(index) = init_mu1;
Mu_2_Vec(index) = init_mu2;
end
plot(Mu_1_Vec)
hold on
plot(Mu_2_Vec)

lung_mean_1=init_mu1;
lung_variance_1=init_segma;

chest_mean_1=init_mu2;
chest_variance_1=segma2_old;

pk_01=init_prob;
pk_11=Prob2_old;

% PLOT GAUSSIAN DISTRIBUTION OF EACH IMAGE AND DISPLAY CHEST AND LUNG IMAGE
% AND THEN CALCULATE THE THRESHOLD THEN COMPARE THE RESULT WITH
 %THE GROUND TRUTH OF THIS IMAGE BY CALCULATING DICE SIMILARITY COEFFECIENT 
%Gaussian distribution 
q=[0:1:255]; 
p_lung_01=pk_01.*(1./sqrt(2.*pi.*lung_variance_1)).*exp (-((q-lung_mean_1).^2)./(2.*lung_variance_1)); 
p_chest_11=pk_11.*(1./sqrt(2.*pi.*chest_variance_1)).*exp (-((q-chest_mean_1).^2)./(2.*chest_variance_1)); 
figure, plot (q,p_lung_01, q,p_chest_11) ,title ( 'Gaussian distribution of first image') 

lung_image_01=zeros(h,w);
chest_image_11=zeros(h,w); 
for kl=1:1:h 
     for k2=1:1:w 
         q=A(kl,k2); 
         if (p_lung_01(q+1)>=p_chest_11(q+1))
             lung_image_01(kl,k2)=q; 
         else
             chest_image_11(kl,k2)=q;
         end
     end
end


imwrite(lung_image_01./255,'lung.bmp','bmp');
imwrite(chest_image_11./255,'chest.bmp','bmp');

%to calculate the dsc for result chest image

Y=double(imread('chest.bmp'));
y1= Y==0 ; %binary image lung and chest in the original 
x1= x>78;         %binary image lung and chest in ground truth 

tp1=y1+x1;         % generation values (0,1,2) 2 mean lung in both            
TP=length(find(tp1==2));% just only to image the common area

%calculating of FN  chest in original  lung in ground truth 
fn=x1-y1;       % (binary Image) to determine the lung in ground truth which isnt in original 
FN=length(find(fn==1));

%calculation of FP  lung in original  , chest in ground truth 
fp=y1-x1;        %binary image to determind the lung in the orignal which isnt in the ground truth 
FP=length(find(fp==1));

DSC=2*TP/(2*TP+FN+FP)



%__________________________________________________________________%

% %get the lung only in the chest image 
% y1_c=(Y<30);    % because the lung is black (zero pixels)
% %get the lung in Ground truth 
% x1_c=(x==255);
% tp_c=x1_c.*y1_c;
% TP_c=length(find(tp_c==1));
% %calculating of FN  chest in original  lung in ground truth 
% fn_c=x1_c-y1_c;
% FN_c=length(find(fn_c==1));
% %calculation of FP  lung in original  , chest in ground truth 
% fp_c=y1_c-x1_c;
% FP_c=length(find(fp_c==1));
% dsc_of_chest_image = 2*TP_c/(2*TP_c+FN_c+FP_c)
% 
% 
% %________________________________________________________________%
% 
% Y=double(imread('lung.bmp'));
% 
% sq=strel('square',3);
% dis=strel('disk',3);
% dis8=strel('disk',8);
% 
% Y=imclose(imopen(Y,dis),dis);
% Y=imclose(imopen(Y,sq),sq);
% Y=imclose(imopen(Y,dis8),dis8);
% 
% %get the lung only in the lung image 
% y1_1=(Y>12);
% 
% %get the lung in Ground truth 
% x1_l=(x==255);
% tp_l=x1_l.*y1_1;
% TP_l=length(find(tp_l==1));
% %calculating of FN  chest in original  lung in ground truth 
% fn_l=x1_l-y1_1;
% FN_l=length(find(fn_l==1));
% %calculation of FP  lung in original  , chest in ground truth 
% fp_l=y1_1-x1_l;
% FP_l=length(find(fp_l==1));
% 
% dsc_of_lung_image = 2*TP_l/(2*TP_l+FN_l+FP_l)
% showing the result (lung , chest ) after testing the model 
% imshow(uint8(lung_image_01)),figure,imshow(uint8(chest_image_11)),figure,plot(q,p_lung_01,q,p_chest_11)
% q1=0:1:255;
% plot(q1,p_lung_01,q1,p_chest_11);
imshow(y1)