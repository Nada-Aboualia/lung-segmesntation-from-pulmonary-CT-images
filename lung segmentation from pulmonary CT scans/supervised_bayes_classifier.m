a1=double(imread('CT_001.bmp'));
a2=double(imread('CT_002.bmp'));
a3=double(imread('CT_003.bmp'));

g1=double(imread('Ground_Truth_CT_001.bmp'));
g2=double(imread('Ground_Truth_CT_002.bmp'));
g3=double(imread('Ground_Truth_CT_003.bmp'));

%reading test image
in_image=double(imread('CT_004.bmp'));
x=imread('Ground_Truth_CT_004.bmp');

A=[a1 a2 a3 ];
G=[g1 g2 g3];


lung =double(((G==255).*A));     %resulting of lung in the origingal image from the grond truth 
chest=double(((G==78).*A));    

%mean and varicance for lung
mu_lung=(sum(nonzeros(lung))/length(nonzeros(lung)));
variance_lung=sum(sum((nonzeros(lung)-mu_lung).^2))/(length(nonzeros(lung))-1);

%mean and variance for chest
mu_chest=(sum(nonzeros(chest))/length(nonzeros(chest)));
variance_chest=sum(sum((nonzeros(chest)-mu_chest).^2))/(length(nonzeros(chest))-1);

pk1=length(nonzeros(lung))/(length(nonzeros(lung))+length(nonzeros(chest))); %perior for lung
pk2=1-pk1;   %perior for chest

q=0:1:255;

%gaussian model for lung
p_lung= pk1.*((1./sqrt(2.*pi.*variance_lung))).*exp(-((q-mu_lung).^2)./(2.*variance_lung));


%gaussian model for chest
p_chest=pk2.*((1./sqrt(2.*pi.*variance_chest))).*exp(-((q-mu_chest).^2)./(2.*variance_chest));

%testing of the model

[r,c]=size(in_image);

lung_image= zeros(r,c);
chest_image= zeros(r,c);

for i=1:r
    for j=1:c
        q =in_image(i,j);
        if(p_lung(q+1))>=(p_chest(q+1))
            lung_image(i,j)=q;
        else
              chest_image(i,j)=q;
        end
    end
end

imwrite(lung_image./255,'lung.bmp','bmp');
imwrite(chest_image./255,'chest.bmp','bmp');

%to calculate the dsc for result chest image

Y=double(imread('chest.bmp'));

%__________________________________________________________________%
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

% showing the result (lung , chest ) after testing the model 
imshow(uint8(lung_image)),figure,imshow(y1),figure,plot(q,p_lung,q,p_chest)
q1=0:1:255;
plot(q1,p_lung,q1,p_chest);



