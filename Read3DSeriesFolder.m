function [img R_img T_img dmeta] = Read3DSeriesFolder(folder)

% Read in images and frame of reference info from a single DICOM series
% stored as single-frame images in a folder.
%
% Requires MATLAB Image Processing toolbox.
%
% Outputs:
%           img  : a 3D pixel array;
%           R_img: a frame of reference in IMAGE SPACE but with PHYSICAL
%                  UNITS (ie. millimetres not pixels) in imref3d format
%           T_img: a 4D transformation matrix suitable for use with imwarp
%                  that maps image space to the DICOM reference coordinate
%                  system
%           dmeta: DICOM metadata, containing only first-level fields whose
%                  value is identical for all instances in the series
%
% This code will assume that there's only one series in the folder, and
% that there are no images with the same Z value. It's therefore not
% appropriate for reading in DWI/DCE myriad other 4D+ data.
%
% Daniel Warren - March 2017

listing = dir(folder);
f = 0;
for j = 1:numel(listing)
    fn = [folder filesep listing(j).name];
    if ~listing(j).isdir && isdicom(fn)
        f = f+1;
        dinfo{f} = dicominfo(fn);
        img(:,:,f) = dicomread(dinfo{f});
    end
end

IOP_xyz = dinfo{1}.ImageOrientationPatient; % direction cosines relating i vector to xyz and j vector to xyz
idir_xyz = IOP_xyz(1:3); % direction vector for rows in xyz
jdir_xyz = IOP_xyz(4:6); % direction vector for columns in xyz
kdir_xyz = cross(idir_xyz,jdir_xyz); % previously compared IPP(1) and IPP(end) but sign problems in sagittal images

idir_xyz = idir_xyz/norm(idir_xyz);
jdir_xyz = jdir_xyz/norm(jdir_xyz);
kdir_xyz = kdir_xyz/norm(kdir_xyz);

img_sl = cellfun(@(x)dot(x.ImagePositionPatient,kdir_xyz),dinfo); % Slice order metric (more reliable than SliceLocation)
img_sl = img_sl-min(img_sl(:));
[img_sl,ord] = sort(img_sl,'ascend');
img = img(:,:,ord);
dinfo = dinfo(ord);

img_vk = cellfun(@(x)x.ImagePositionPatient,dinfo,'UniformOutput',false); % img_vk will be vector of pixel centres in k
img_vk = bsxfun(@minus,[img_vk{:}],dinfo{1}.ImagePositionPatient);
img_vk = sqrt(sum(img_vk.^2,1));

img_di = dinfo{1}.PixelSpacing(1); % DICOM specification says first element of PS is rows
img_dj = dinfo{1}.PixelSpacing(2); % and second element is columns
img_dk = img_vk(2)-img_vk(1);

img_vi = img_di*double(0:(-1+dinfo{1}.Rows)); % img_vi is vector of pixel centres in i
img_vj = img_dj*double(0:(-1+dinfo{1}.Columns)); % img_vj is vector of pixel centres in j

IPP_xyz = dinfo{1}.ImagePositionPatient; % (x,y,z) position of (i,j)=(0,0)

rotmat = [idir_xyz jdir_xyz kdir_xyz]; % Rotation matrix to map ijk directions to patient xyz directions
rotmat(4,4) = 1;
tmat = [[eye(3); IPP_xyz'] [0;0;0;1]]; % Translation matrix to map origin to IPP(0)

img_fovi = [img_vi(1)-img_di/2 img_vi(end)+img_di/2];
img_fovj = [img_vj(1)-img_dj/2 img_vj(end)+img_dj/2];
img_fovk = [img_vk(1)-img_dk/2 img_vk(end)+img_dk/2];

R_img = imref3d(size(img),img_fovj,img_fovi,img_fovk); % Frame of reference for image space
T_img = rotmat'*tmat; % Transformation matrix from image coordinates to real space
                     % Matrix is right-multiplied by imwarp, so earlier
                     % operations should appear first

% Keep a subset of the DICOM metadata - only first-level fields with numeric
% or char datatypes, whose value is the same for every instance
dmeta = struct();
fn = fieldnames(dinfo{1});
for i = 1:numel(fn)
    keep = false;
    refinfo = dinfo{1}.(fn{i});
    for j = 2:numel(dinfo)
        if isfield(dinfo{j},fn{i})
            testinfo = dinfo{j}.(fn{i});
            if strcmp(class(refinfo),class(testinfo)) && ... 
                ( ...
                 ( isnumeric(refinfo) && (numel(refinfo) == numel(testinfo)) && all(refinfo == testinfo) ) ...
                || ...
                 ( ischar(refinfo) && strcmp(refinfo,testinfo) ) ...
                )
                    keep = true;
            end
        else
            keep = false;
        end
        if ~keep; break; end;
    end
    if keep
        dmeta.(fn{i}) = dinfo{1}.(fn{i});
    end
end

end