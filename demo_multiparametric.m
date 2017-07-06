% This demo script downloads some sample multiparametric MRI data,
% calculates apparent diffusion coefficient (ADC) from a diffusion-weighted
% imaging (DWI) sequence, and then displays a GUI where fused T1, T1+Gd
% enhancement and ADC can be visualised for a chosen slice.
%
% Daniel R Warren
% June 2017
% http://github.com/drw25

%% Download and extract some sample data

disp('Downloading data...');

% Images from giveascan.org, patient P40002, study 2011-03-14
url_t1 = 'http://www.giveascan.org/download/?items=12305'; % T1-weighted
url_t1c = 'http://www.giveascan.org/download/?items=12300'; % T1-weighted + contrast
url_dwi = 'http://www.giveascan.org/download/?items=12299'; % DWI

fn_t1 = [tempdir filesep 'multimodal_t1.zip'];
fn_t1c = [tempdir filesep 'multimodal_t1c.zip'];
fn_dwi = [tempdir filesep 'multimodal_dwi.zip'];

urlwrite(url_t1,fn_t1);
urlwrite(url_t1c,fn_t1c);
urlwrite(url_dwi,fn_dwi);

folder_t1 = [tempdir filesep 'multimodal_t1'];
folder_t1c = [tempdir filesep 'multimodal_t1c'];
folder_dwi = [tempdir filesep 'multimodal_dwi'];

mkdir(folder_t1);
mkdir(folder_t1c);
mkdir(folder_dwi);

unzip(fn_t1,folder_t1);
unzip(fn_t1c,folder_t1c);
unzip(fn_dwi,folder_dwi);

delete(fn_t1);
delete(fn_t1c);
delete(fn_dwi);

%% Read in images

disp('Reading sample data...');

[img_prim R_prim T_prim] = Read3DSeriesFolder(folder_t1);
[img_t1c R_t1c T_t1c] = Read3DSeriesFolder(folder_t1c);

% in the DWI dataset, Siemens b-value fields have been blanked out by the
% anonymizer, but fortunately values are given as part of the SequenceName
[img_dwi bvals R_dwi T_dwi] = Read4DSeriesFolder(folder_dwi,'SequenceName');
bvals = cellfun(@(x)sscanf(x,'*ep_b%d'),bvals);

% cast data to double to facilitate calculations
img_prim = double(img_prim);
img_t1c = double(img_t1c);
img_dwi = double(img_dwi);

%% Delete sample data files

delete([folder_t1 filesep '*']);
rmdir(folder_t1);
delete([folder_t1c filesep '*']);
rmdir(folder_t1c);
delete([folder_dwi filesep '*']);
rmdir(folder_dwi);

%% Calculate ADC from DWI

disp('Calculating ADC...');

% Assumed relation for ADC is S = S0*exp(-ADC*bval)

bvals = reshape(bvals,[1 1 1 numel(bvals)]);

S0 = img_dwi(:,:,:,bvals == 0);
logratio = -log(bsxfun(@rdivide,img_dwi,S0));

% find slope of logratio vs bvals by simple linear regression
SXY = sum(bsxfun(@times,bvals,logratio),4)-numel(bvals)*mean(bvals,4)*mean(logratio,4);
SXX = sum(bvals.^2,4)-numel(bvals)*mean(bvals,4).^2;

img_adc = SXY/SXX;
img_adc(~isfinite(img_adc)) = 0;

clear logratio S0 SXY SXX;

%% Map images on to T1

disp('Transforming images...');

T_c_to_prim = T_t1c/T_prim;
T_c_to_prim(:,4) = [0;0;0;1]; % fix precision errors
tform_c_to_prim = affine3d(T_c_to_prim);
img_c_to_prim = imwarp(img_t1c,R_t1c,tform_c_to_prim,'OutputView',R_prim);

T_dwi_to_prim = T_dwi/T_prim;
T_dwi_to_prim(:,4) = [0;0;0;1]; % fix precision errors
tform_dwi_to_prim = affine3d(T_dwi_to_prim);
img_adc_to_prim = imwarp(img_adc,R_dwi,tform_dwi_to_prim,'OutputView',R_prim);
img_b0_to_prim = imwarp(img_dwi(:,:,:,bvals==0),R_dwi,tform_dwi_to_prim,'OutputView',R_prim);

%% Calculate 3D colour maps for each sequence
% (this allows simpler callback functions in the GUI cf. demo_contours)

disp('Calculating colour maps..,')

constrain = @(x) x.*(x > 0 & x <= 1) + 1*(x > 1); % cap x in range 0-1
applycmap = @(x,c) reshape(c(1+round((size(c,1)-1)*x),:),[size(x) 3]); % make colour mapped image from 2D array with values in range 0-1
ice = @(n) flipud(1-hot(n)); % cold colour map

t1_cols = gray(256);
t1_scaled = constrain(img_prim/750);
t1_map = applycmap(t1_scaled,t1_cols);

adc_cols = ice(256);
adc_scaled = constrain((img_adc_to_prim-1e-3)/2e-3);
adc_map = applycmap(adc_scaled,adc_cols);

t1en_cols = hot(256);
t1en_scaled = constrain((img_c_to_prim-img_prim)/750);
t1en_scaled(~isfinite(t1en_scaled)) = 0;
t1en_map = applycmap(t1en_scaled,t1en_cols);

% alpha maps - use pixel max colour value as default
adc_alpha = max(adc_map,[],4);
t1en_alpha = max(t1en_map,[],4);

%% Create GUI

disp('Creating GUI...');

f = figure('Units','Pixels','Position',[100 100 600 672],'Name','Move sliders to change slice/alpha');
f_lab = uicontrol('Units','Pixels','Position',[44 642 512 30],'Style','text','String','Background = T1; Hot = T1+Gd enhancement; Cold = ADC from DWI');
a = axes('Units','Pixels','Position',[44 130 512 512]);
ci = uicontrol('Units','Pixels','Position',[88 22 468 36],'Style','slider','Value',0.8);
ci_lab = uicontrol('Units','Pixels','Position',[22 22 66 36],'Style','text','String','Slice');
cm1 = uicontrol('Units','Pixels','Position',[88 58 468 36],'Style','slider','Value',0.5);
cm1_lab = uicontrol('Units','Pixels','Position',[22 58 66 36],'Style','text','String','ADC alpha');
cm2 = uicontrol('Units','Pixels','Position',[88 94 468 36],'Style','slider','Value',0.5);
cm2_lab = uicontrol('Units','Pixels','Position',[22 94 66 36],'Style','text','String','T1 enh alpha');

% alpha blend for opaque background (https://en.wikipedia.org/wiki/Alpha_compositing#Alpha_blending) 
alphablend = @(rgb1,rgb2,a1,a2) bsxfun(@times,rgb2,a2)+bsxfun(@times,rgb1,1-a2);

% slice color map array
sl = @(x,i) squeeze(x(:,:,i,:));

% chained functions for callback
blendslice1 = @(i,mod1) alphablend(sl(t1_map,i),sl(adc_map,i),1,mod1*sl(adc_alpha,i));
blendslice2 = @(i,mod1,mod2) alphablend(blendslice1(i,mod1),sl(t1en_map,i),1,mod2*sl(t1en_alpha,i));

sz = size(t1_map);
callback = @(varargin) imshow(blendslice2(round(1+get(ci,'Value')*(sz(3)-1)),get(cm1,'Value'),get(cm2,'Value')));

set(ci,'Callback',callback);
set(cm1,'Callback',callback);
set(cm2,'Callback',callback);

callback();