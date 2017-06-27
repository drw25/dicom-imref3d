% This demo script downloads two sample MRI datasets with different
% orientations, reads them in using Read3DSeriesFolder, and then fuses them
% on a common coordinate system.
%
% Daniel Warren - June 2017

%% Download some sample data (server is a bit slow)

disp('Downloading data...');

dir1 = [tempdir filesep 'Series1'];
dir2 = [tempdir filesep 'Series2'];

mkdir(dir1);
mkdir(dir2);

f = ftp('medical.nema.org');
cd(f,'medical/dicom/DataSets/WG16/Philips/ClassicSingleFrame/Knee/DICOM');

h = waitbar(0,'Downloading series 1...');
for i = 27:50
    mget(f,['IM_' num2str(i,'%04d')],dir1);
    waitbar((i-26)/24,h);
end
close(h);

h = waitbar(0,'Downloading series 2...');
for i = 53:76
    mget(f,['IM_' num2str(i,'%04d')],dir2);
    waitbar((i-52)/24,h);
end
close(h);

close(f);

%% Read in images

disp('Reading series...');

[img1 R1 T1 meta1] = Read3DSeriesFolder(dir1);
[img2 R2 T2 meta2] = Read3DSeriesFolder(dir2);

%% Delete sample data files

delete([dir1 filesep '*']);
rmdir(dir1);
delete([dir2 filesep '*']);
rmdir(dir2);

%% Map image 2 on to image 1

disp('Transforming image...');

T2to1 = T2/T1;
T2to1(:,4) = [0;0;0;1]; % fix precision errors
tform2to1 = affine3d(T2to1);
img2to1 = imwarp(img2,R2,tform2to1,'OutputView',R1);

%% Display output

figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);

subplot(1,3,1);
imshow(img1(:,:,round(end/2)),[0 500],'ColorMap',gray(256));
title(['Series' num2str(meta1.SeriesNumber)]);

subplot(1,3,2);
imshow(img2(:,:,round(end/2)),[0 500],'ColorMap',gray(256));
title(['Series' num2str(meta2.SeriesNumber)]);

subplot(1,3,3);
imshowpair(img1(:,:,round(end/2)),img2to1(:,:,round(end/2)));
title('Fused');

disp('Done.');