function [ mask ] = Mask_Circle(D,s )
%MASK_CIRCLE Creates a circular mask of diameter D and dimensions [s,s]
%[ mask ] = Mask_Circle(D,s )
% D     -   Diameter of the mask
% s     -   Dimension of the mask array
% mask  -   Binary array with dimensions [s,s]. mask=1 inside a circle of
%           diameter D circumscript in the mask array
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

% Define default mask diameter
if nargin<2||isempty(s),s=[D D];end
if size(s,2)==1,s=[s s];end

% Create arrays of horizontal and vertical coordinates of the image
[X,Y]=meshgrid(1:s(1),1:s(2));
cx=ceil(s(1)/2);
cy=ceil(s(2)/2);

%Define an array with (radial) distances to the pixel with coordinates
% (cx,cy);
circle=sqrt((X-cx).^2+(Y-cy).^2);

% Define circular mask. mask=1 for all points with radial distance to
% (cx,cy) equal or less to the circlar mask radius
mask=circle<=(D/2);

end

