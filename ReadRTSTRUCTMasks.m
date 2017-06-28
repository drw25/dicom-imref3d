function [masks names img R_img T_img] = ReadRTSTRUCTMasks(file,folder)

% Read structure masks from a DICOM RTSTRUCT file, where contours have been
% defined on a DICOM image series stored in a folder.
%
% Requires MATLAB Image Processing toolbox.
%
% Outputs:
%           masks: a 4D logical array labelling whether each voxel is
%                  inside each structure
%           names: structure names (same order as dimension 4
%                              of masks)
%           img  : image data, as Read3DSeriesFolder
%           R_img: image and mask frame of reference, as Read3DSeriesFolder
%           T_img: image and mask transformations, as Read3DSeriesFolder
%
% Folder should contain one simple image series, as for Read3DSeriesFolder.
%
% Daniel Warren - June 2017

[img R_img T_img dmeta] = Read3DSeriesFolder(folder); % Image series

dinfo = dicominfo(file); % RTSTRUCT contents

% Extract ROI numbers and names for each structure
fn = fieldnames(dinfo.StructureSetROISequence);
names = cell(1,numel(fn));
roinum = zeros(1,numel(fn));
for i = 1:numel(fn)
    item = dinfo.StructureSetROISequence.(fn{i});
    names{i} = item.ROIName;
    roinum(i) = item.ROINumber;
end

% Extract ROI coordinate data for each structure
contourdata = cell(1,numel(names));
fn = fieldnames(dinfo.ROIContourSequence);
for i = 1:numel(fn)
    item = dinfo.ROIContourSequence.(fn{i});
    refroi = item.ReferencedROINumber;
    if isfield(item,'ContourSequence')
        contourdata{refroi == roinum} = item.ContourSequence;
    end
end

% For each ROI segment, find the appropriate image slice, convert the segment
% coordinates to the image coordinate system and rasterize the contour.
% Store the result in the 4D variable masks.

masks = false([size(img) numel(names)]);
for i = 1:numel(names)
    if ~isempty(contourdata{i})
        fn = fieldnames(contourdata{i});
        for j = 1:numel(fn)

            item = contourdata{i}.(fn{j});
            refuid = item.ContourImageSequence.Item_1.ReferencedSOPInstanceUID;

            sliceind = find(strcmp(refuid,dmeta.SOPInstanceUID),1,'first');

            if isempty(sliceind)
                warning(['Missing image with UID reference in RTSTRUCT: ' refuid]);
            end

            contourpts = item.ContourData;

            % Find conversion from RCS (used in RTSTRUCT) to DICOM pixel indices by
            % inverting the transform given in Equation C.7.6.2.1-1 of the DICOM
            % specification (including an extra 1 because otherwise the matrix is
            % singular - this is OK because 3rd component of input\output vector is zero)

            IOP = dmeta.ImageOrientationPatient{sliceind};
            IPP = dmeta.ImagePositionPatient{sliceind};
            PS = dmeta.PixelSpacing{sliceind};

            M = [IOP(1)*PS(1) IOP(4)*PS(2) 0 IPP(1);
                 IOP(2)*PS(1) IOP(5)*PS(2) 0 IPP(2);
                 IOP(3)*PS(1) IOP(6)*PS(2) 1 IPP(3);
                 0            0            0 1      ];

            contourpts = reshape(contourpts,[3 size(contourpts,1)/3]);
            contourpts(4,:) = 1;

            contourpts_rcs = zeros(size(contourpts));
            for k = 1:size(contourpts,2)
                contourpts_rcs(:,k) = M\contourpts(:,k);
            end

            xpoly = 1+contourpts_rcs(1,:); % since DICOM index is from 0, but MATLAB
            ypoly = 1+contourpts_rcs(2,:); % works from 1

            contourmask = poly2mask(xpoly,ypoly,size(img,1),size(img,2));
            masks(:,:,sliceind,i) = contourmask | masks(:,:,sliceind,i); % OR in case multiple segments on each slice
        end     
    end
end

end
