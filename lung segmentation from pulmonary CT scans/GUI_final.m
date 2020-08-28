
function varargout = GUI_final(varargin)
% GUI_FINAL MATLAB code for GUI_final.fig
%      GUI_FINAL, by itself, creates a new GUI_FINAL or raises the existing
%      singleton*.
%
%      H = GUI_FINAL returns the handle to a new GUI_FINAL or the handle to
%      the existing singleton*.
%
%      GUI_FINAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FINAL.M with the given input arguments.
%
%      GUI_FINAL('Property','Value',...) creates a new GUI_FINAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_final_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_final_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_final

% Last Modified by GUIDE v2.5 14-Aug-2020 07:50:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_final_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_final_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_final is made visible.
function GUI_final_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_final (see VARARGIN)

% Choose default command line output for GUI_final
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_final wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_final_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%--- %Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Filename1 ] = uigetfile('*.bmp','File Selector');
[Filename2 ] = uigetfile('*.bmp','File Selector');
[Filename3 ] = uigetfile('*.bmp','File Selector');
OR1=strcat( Filename1) ;
OR2=strcat( Filename2) ;
OR3=strcat( Filename3) ;
OR= [OR1 OR2 OR3];

a1=imread(strcat (OR1));
a2=imread(strcat( OR2));
a3=imread(strcat( OR3));
A1 = [a1 a2 a3];
axes(handles.axes1);
imshow(A1);
set(handles.edit1,'string',OR);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Filename4] = uigetfile('*.bmp','File Selector');
OR4=strcat( Filename4)
a4= imread(OR4);
axes(handles.axes5);
imshow(a4);
set(handles.edit4,'string',OR4);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Filename5 ] = uigetfile('*.bmp','File Selector');
[Filename6 ] = uigetfile('*.bmp','File Selector');
[Filename7 ] = uigetfile('*.bmp','File Selector');
GT1=strcat( Filename5) ;
GT2=strcat( Filename6) ;
GT3=strcat( Filename7) ;
GT= [GT1 GT2 GT3];

g1=imread(strcat (GT1));
g2=imread(strcat( GT2));
g3=imread(strcat( GT3));
G1 = [g1 g2 g3];
axes(handles.axes2);
imshow(G1);
set(handles.edit2,'string',GT);
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Filename8] = uigetfile('*.bmp','File Selector');
GT4=strcat(Filename8);
g4= imread(GT4);
axes(handles.axes6);
imshow(g4);
set(handles.edit3,'string',GT4);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
n=get(handles.popupmenu1,'value');
if n==2
Y= double([getimage(handles.axes1) getimage(handles.axes5)]);
x= [getimage(handles.axes2) getimage(handles.axes6)];

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

%calculating TP   lung in both (origninal and ground truth)
tp1=y2+x1;         % generation values (0,1,2) 2 mean lung in both 
tp=tp1>1;          % just only to image the common area 
TP=length(find(tp==1));

%calculating of FN  chest in original  lung in ground truth 
fn=x1-y2;       % (binary Image) to determine the lung in ground truth which isnt in original 
FN=length(find(fn==1));

%calculation of FP  lung in original  , chest in ground truth 
fp=y2-x1;        %binary image to determind the lung in the orignal which isnt in the ground truth 
FP=length(find(fp==1));

Dice =2*TP/(2*TP+FN+FP)
set(handles.edit5,'string',Dice)
axes(handles.axes3)
imshow(y2)


elseif n==3
%reading test image

A = double(getimage(handles.axes1));
G = double(getimage(handles.axes2));
in_image= double(getimage(handles.axes5));
x= getimage(handles.axes6);

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

Dice=2*TP/(2*TP+FN+FP)

% showing the result (lung , chest ) after testing the model 
axes(handles.axes3);
imshow(y1);
set(handles.edit5,'string',num2str(Dice));
axes(handles.axes4);
plot(q,p_lung,q,p_chest)
q1=0:1:255;
plot(q1,p_lung,q1,p_chest);

elseif n==4 

init_mu1 = 30;
init_segma = 100;
init_prob = 0.5;
init_mu2 = 200;
segma2_old = 100;
Prob2_old = 0.5;
 
A= double([getimage(handles.axes1) getimage(handles.axes5)]);
x= ([getimage(handles.axes2) getimage(handles.axes6)]);
[h w] = size(A);


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
% axes(handles.axes4)
% plot(Mu_1_Vec)
% hold on
% plot(Mu_2_Vec)

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
axes(handles.axes4);
plot (q,p_lung_01, q,p_chest_11);
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
%showing results
Dice =2*TP/(2*TP+FN+FP);
set(handles.edit5,'string',num2str(Dice));
axes(handles.axes3);
imshow(y1);
   
end


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes1,'reset');
cla(handles.axes2,'reset');
cla(handles.axes3,'reset');
cla(handles.axes4,'reset');
cla(handles.axes5,'reset');
cla(handles.axes6,'reset');
set(handles.edit1,'String','Training img names');
set(handles.edit2,'String','GT of trainig img names');
set(handles.edit3,'String','Test image name');
set(handles.edit4,'String','GT of test img name');
set(handles.edit5,'String','Dice similariry');
