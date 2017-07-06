% This demo script downloads a sample DICOM-RT dataset, reads in the CT and
% structure set using ReadRTSTRUCTMasks, and then displays a GUI where
% structure contours can be visualised on top of the CT for a chosen slice.
%
% Daniel R Warren
% June 2017
% http://github.com/drw25

%% Download and extract some sample data

disp('Downloading data...');

url = ['http://slicer.kitware.com/midas3/download/item/10704/EclipseEntDicomRt.zip'];

fn_zip = [tempdir filesep 'rtstruct.zip'];

urlwrite(url,fn_zip);

folder = [tempdir filesep 'rtstruct'];
mkdir(folder);

unzip(fn_zip,folder);
delete(fn_zip);

imgfolder = [folder filesep 'img'];
mkdir(imgfolder);
movefile([folder filesep 'CT*'],imgfolder);

listing = dir(folder);
structfile = listing(~cellfun(@isempty,strfind({listing.name},'RS.'))).name;

%% Read image and masks

disp('Reading sample data...');

[masks names img R_img T_img] = ReadRTSTRUCTMasks([folder filesep structfile],imgfolder);

%% Delete sample data files

delete([imgfolder filesep '*']);
rmdir(imgfolder);
delete([folder filesep '*']);
rmdir(folder);

%% Create GUI

disp('Creating GUI...');

f = figure('Units','Pixels','Position',[100 100 600 600],'Name','Move slider to change slice');
a = axes('Units','Pixels','Position',[44 58 512 512]);
c = uicontrol('Units','Pixels','Position',[44 22 512 36],'Style','slider');

sz = size(masks);

cols = [0 0 0; hsv(sz(4))]; % colours for structure contours (first colour must always be black for compositeslice to work)
cols(1:2:end,:) = (2/3)*cols(1:2:end,:);

stack = @(x) cat(3,x{:});

% Chain of functions for callback
masks2cell = @(i) squeeze(mat2cell(masks(:,:,i,:),sz(1),sz(2),1,ones(1,sz(4)))); % masks on slice i as cell
bwcontours2cell = @(i) cellfun(@bwperim,masks2cell(i),'UniformOutput',false); % contours on slice i as cell
topcontour = @(i) max(bsxfun(@times,stack(bwcontours2cell(i)),permute(1:sz(4),[3 1 2])),[],3); % greatest contour index on slice i as array
contourslice = @(i) reshape(cols(1+topcontour(i),:),[sz(1) sz(2) 3]); % colour map showing greatest index contour
imgslice = @(i) repmat(0.5+(-1000+double(img(:,:,i)))/400,[1 1 3]); % colour map of image slice, hardcoded window [-200; 200], assuming rescale intercept of -1000
compositeslice = @(i) bsxfun(@times,~any(stack(bwcontours2cell(i)),3),imgslice(i))+contourslice(i); % colour map composite of image and contour

showslice = @(i) imshow(compositeslice(i));
printlabels = @(ax) cellfun(@(x,y,s,c) text(ax,x,y,s,'Units','pixels','FontSize',8,'Color',c), ... % print coloured list of structure names
                    num2cell(12*ones(1,sz(4))), num2cell(512-12*(1:sz(4))), names, ...
                    mat2cell(cols(2:end,:),ones(1,sz(4)),3)');

showlabelslice = @(i) printlabels(get(showslice(i),'Parent')); % print labels on composite slice image

set(c,'CallBack',@(h,e) showlabelslice(round(1+get(h,'Value')*(sz(3)-1))));
set(c,'Value',0.75);
showlabelslice(round(1+0.75*(sz(3)-1)));

disp('Done.');