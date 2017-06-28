function [img d4vals R_img T_img dmeta] = Read4DSeriesFolder(folder,d4tag)

% Read in 4D images and frame of reference info from a single DICOM series
% stored as single-frame images in a folder.
%
% Values for dimension 4 are taken from the DICOM tag specified by d4tag,
% which should be a char containing the name of a field from the output
% structure of dicominfo.
%
% Examples of previously useful d4tag values are:
%  'FrameReferenceTime' : frame time in GE dynamic PET
%  'RepetitionTime': input to T1 map calculations in MRI
%  'TriggerTime': frame time in GE DCE-MRI
%  'Private_0019_1024' : frame time in GE perfusion CT
%  'Private_0019_100c' : b-value in Siemens diffusion-weighted MRI
%  'Private_0043_1039' : b-value in GE diffusion-weighted MRI
%
% Code is essentially the same as Read3DSeriesFolder, except with an
% additional sorting step for the 4th dimension variable.
%
% Requires MATLAB Image Processing toolbox.
%
% Outputs:
%           img   : a 4D pixel array
%           d4vals: a vector of values for the 4th dimension variable. Note
%                   that some numerical values may not be in a numerical
%                   datatype eg. if the DICOM tag is not in MATLAB's
%                   dictionary and the VR is IS/DS
%           R_img : a frame of reference in IMAGE SPACE but with PHYSICAL
%                   UNITS (ie. millimetres not pixels) in imref3d format
%           T_img : a 4D transformation matrix suitable for use with imwarp
%                   that maps image space to the DICOM reference coordinate
%                   system
%           dmeta : Subset of DICOM metadata. UIDs and 3 important geometric
%                   fields are retained for each slice (IOP, IPP, PS).
%                   Remainder is only first-level fields whose value is
%                   identical for all instances in the series
%
% This code will assume that there's only one series in the folder.
%
% Daniel Warren - June 2017

% Initially read the data+metadata into 3D+1D  '_lin' variables
listing = dir(folder);
f = 0;
for j = 1:numel(listing)
    fn = [folder filesep listing(j).name];
    if ~listing(j).isdir && isdicom(fn)
        f = f+1;
        dinfo_lin{f} = dicominfo(fn); % metadata (cell since struct fields can vary)
        img_lin(:,:,f) = dicomread(dinfo_lin{f}); % pixeldata
    end
end

IOP_xyz = dinfo_lin{1}.ImageOrientationPatient; % direction cosines relating i vector to xyz and j vector to xyz
idir_xyz = IOP_xyz(1:3); % direction vector for rows in xyz
jdir_xyz = IOP_xyz(4:6); % direction vector for columns in xyz
kdir_xyz = cross(idir_xyz,jdir_xyz); % previously compared IPP(1) and IPP(end) but sign problems in sagittal images

idir_xyz = idir_xyz/norm(idir_xyz);
jdir_xyz = jdir_xyz/norm(jdir_xyz);
kdir_xyz = kdir_xyz/norm(kdir_xyz);

img_sl_lin = cellfun(@(x)dot(x.ImagePositionPatient,kdir_xyz),dinfo_lin); % Slice order metric (more reliable than SliceLocation)
img_sl_lin = img_sl_lin-min(img_sl_lin(:));

% Find values for 4th dimension variable
d4vals_lin = cellfun(@(x)x.(d4tag),dinfo_lin,'UniformOutput',false);
if all(cellfun(@isnumeric,d4vals_lin)) && (all(cellfun(@numel,d4vals_lin) == 1))
    d4vals_lin = [d4vals_lin{:}];
    d4num = true;
elseif isa(d4vals_lin{1},'uint8')
    d4vals_lin = cellfun(@(x)char(x'),d4vals_lin,'UniformOutput',false); % assume that uint8 datatypes are actually char
    d4num = false;
else
    d4num = false;
end

img_sl = sort(unique(img_sl_lin),'ascend');
d4vals = sort(unique(d4vals_lin));

% Reshape data+metadata into 4D+2D variables
img = zeros([size(img_lin,1) size(img_lin,2) numel(img_sl) numel(d4vals)],class(img_lin));
dinfo = cell([numel(img_sl) numel(d4vals)]);
for i = 1:numel(img_sl)
    for j = 1:numel(d4vals)
        if d4num
            ind = img_sl_lin == img_sl(i) & d4vals_lin == d4vals(j);
        else
            ind = img_sl_lin == img_sl(i) & strcmp(d4vals_lin,d4vals{j});
        end
        img(:,:,i,j) = img_lin(:,:,ind);
        dinfo{i,j} = dinfo_lin{ind};
    end
end

clear img_lin dinfo_lin d4vals_lin;

% Calculate image frame of reference

img_vk = cellfun(@(x)x.ImagePositionPatient,dinfo(:,1),'UniformOutput',false); % img_vk will be vector of pixel centres in k
img_vk = bsxfun(@minus,[img_vk{:}],dinfo{1,1}.ImagePositionPatient);
img_vk = sqrt(sum(img_vk.^2,1));

img_di = dinfo{1,1}.PixelSpacing(1); % DICOM specification says first element of PS is rows
img_dj = dinfo{1,1}.PixelSpacing(2); % and second element is columns
img_dk = img_vk(2)-img_vk(1);

img_vi = img_di*double(0:(-1+dinfo{1,1}.Rows)); % img_vi is vector of pixel centres in i
img_vj = img_dj*double(0:(-1+dinfo{1,1}.Columns)); % img_vj is vector of pixel centres in j

IPP_xyz = dinfo{1,1}.ImagePositionPatient; % (x,y,z) position of (i,j)=(0,0)

rotmat = [idir_xyz jdir_xyz kdir_xyz]; % Rotation matrix to map ijk directions to patient xyz directions
rotmat(4,4) = 1;
tmat = [[eye(3); IPP_xyz'] [0;0;0;1]]; % Translation matrix to map origin to IPP(0)

img_fovi = [img_vi(1)-img_di/2 img_vi(end)+img_di/2];
img_fovj = [img_vj(1)-img_dj/2 img_vj(end)+img_dj/2];
img_fovk = [img_vk(1)-img_dk/2 img_vk(end)+img_dk/2];

R_img = imref3d([size(img,1) size(img,2) size(img,3)],img_fovj,img_fovi,img_fovk); % Frame of reference for image space
T_img = rotmat'*tmat; % Transformation matrix from image coordinates to real space
                     % Matrix is right-multiplied by imwarp, so earlier
                     % operations should appear first

dmeta = struct();
                     
% Keep a subset of DICOM metadata - only first-level fields with numeric or
% char datatypes, whose value is the same for every instance
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

% Keep some important DICOM metadata for every slice

fieldcontents = @(f)cellfun(@(x){x.(f)},dinfo);

dmeta.SOPInstanceUID = fieldcontents('SOPInstanceUID');
dmeta.ImageOrientationPatient =  fieldcontents('ImageOrientationPatient');
dmeta.ImagePositionPatient = fieldcontents('ImagePositionPatient');
dmeta.PixelSpacing = fieldcontents('PixelSpacing');

end