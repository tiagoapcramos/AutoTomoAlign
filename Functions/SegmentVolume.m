function [ volume_segmented ] = SegmentVolume(volume,method )
% SEGMENTVOLUME Segments reconstructed volume into a binary array 
%[ volume_segmented ] = SegmentVolume(volume,method )
%
% The binary segmentation is performed using the Otsu's method
% volume           -  reconstructed volume. Sample voxels assumed positive
%                     (For real data applications, phase-contrast
%                     projections can be negative. For those cases add a
%                     minus sign '-' to the input volume.
% method           -  'slice' for a slice by slice segmentation 'volume'
%                     for segmentation of the whole volume at once
% volume_segmented - Segmented volume: 1 for object, 0 for air/background
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

% Define default segmentation method as a slice-by-slice approach
if nargin<2||isempty(method),method='slice';end

if strcmp(method,'volume')
    % Compute a single threshold level on volume (see Otsu's method)
    t=multithresh(volume,1);
    % Binarizes/Segments volume according to intensity threshold t
    volume_segmented=imquantize(volume,t);
else
    volume_segmented=volume;
    for k=1:size(volume_segmented,3);
        % Compute a single threshold level on a slice from volume
        t=multithresh(volume_segmented(:,:,k),1);
        % Binarizes/Segments slice according to intensity threshold t
        volume_segmented(:,:,k)=imquantize(volume_segmented(:,:,k),t);
    end
end

% Adjust intensity levels: 1 for sample 0 for air
volume_segmented=volume_segmented-1;

end

