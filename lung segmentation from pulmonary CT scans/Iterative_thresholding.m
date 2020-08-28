Y=double(imread('CT_002.bmp'));
x=imread('Ground_Truth_CT_002.bmp');
T= 128; % Initial Threshold
for i=1:10
Mu1=mean(nonzeros(((Y>T).*Y))); % Compute the mean of the nonzero values
% in the matrix [(Y>T).*Y]
Mu2=mean(nonzeros(((Y<=T).*Y))); % Compute the mean of the nonzero values
% in the matrix [(Y<=T)*Y]
T_updated= (Mu1 + Mu2)/2 % Updated value of threshold
if T_updated==T
break
else
T=T_updated;
end
end
y2 = Y<T_updated; %binary image lung and chest in the original 
x1= x>78;         %binary image lung and chest in ground truth 

sq=strel('square',3);
dis=strel('disk',3);
dis8=strel('disk',8);

y1=imclose(imopen(y2,dis),dis);
y1=imclose(imopen(y1,sq),sq);
y1=imclose(imopen(y1,dis8),dis8);



%calculating TP   lung in both (origninal and ground truth)
tp1=y1+x1;         % generation values (0,1,2) 2 mean lung in both 
tp=tp1>1;          % just only to image the common area 
TP=length(find(tp==1));

%calculating of FN  chest in original  lung in ground truth 
fn=x1-y1;       % (binary Image) to determine the lung in ground truth which isnt in original 
FN=length(find(fn==1));

%calculation of FP  lung in original  , chest in ground truth 
fp=y1-x1;        %binary image to determind the lung in the orignal which isnt in the ground truth 
FP=length(find(fp==1));

DSC=2*TP/(2*TP+FN+FP)


% to illustrate the diference in images 
% imshow(x1),figure,imshow(y1),figure,imshow(tp)
% imshow(x1),figure,imshow(y1),figure,imshow(fn)
% imshow(x1),figure,imshow(y1),figure,imshow(fp)