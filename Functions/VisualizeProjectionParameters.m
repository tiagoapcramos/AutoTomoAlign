function [ h ] = VisualizeProjectionParameters(proj_geom )
%VISUALIZEPROJECTIONPARAMETERS Shows projection geometry setup
% [ h ] = VisualizeProjectionParameters(proj_geom )
%Creates a 3D visual representation of the projections orientation
%relative to the "sample". The sample is represented by a cube centered in
%the origin, and the projection images by squares.
%
%
% proj_geom     - Structure array with 3D parallel beam geometry (see
%                 ASTRA's toolbox documentation)
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

%Define a default distance from the origin to the projections centre
distance=8;

%Define number of projections
num_proj=size(proj_geom.Vectors,1);

%Determine sizes of the sample (All geometries are normalized to the sample
%dimensions. The centre of the detector pixel needs to be normalized
N=proj_geom.DetectorColCount;
M=proj_geom.DetectorRowCount;

%Allocate variable to define projection planes position
%(Array dimensions -> [4 corners, 3 coordinates, n_angles planes])
planes=zeros(4,3,num_proj);

%Define detector centre coordinates + offset + position of corners
%Loop over different projections
for i=1:num_proj
    
    %Define centre of the projection plane according to angular
    %orientation. (centre + offset=distance)
    planes(:,:,i)=repmat(proj_geom.Vectors(i,4:6)./[N,N,M]+...
        distance*proj_geom.Vectors(i,1:3),[4,1]); 
    u=proj_geom.Vectors(i,7:9);
    v=proj_geom.Vectors(i,10:12);
    
    %Define coordinates of the projection corners from its central point
    planes(:,:,i)=planes(:,:,i)+0.5*[u+v;-u+v;-u-v;u-v];
end

%Create new figure
h=figure('color','white'); hold on

%Define coordinates of the cube corners
x=[-1 1 1 -1 -1 -1;1 1 -1 -1 1 1;1 1 -1 -1 1 1;-1 1 1 -1 -1 -1]*.5;
y=[-1 -1 1 1 -1 -1;-1 1 1 -1 -1 -1;-1 1 1 -1 1 1;-1 -1 1 1 1 1]*.5;
z=[-1 -1 -1 -1 -1 1;-1 -1 -1 -1 -1 1;1 1 1 1 -1 1;1 1 1 1 -1 1]*.5;

%Images the 3D central cube (sample)
for i=1:6
    h=patch(x(:,i),y(:,i),z(:,i),'k');
    set(h,'edgecolor','w')
end

%Define an array with different RGB colours for the projection planes
color=hsv(num_proj);

%Loop over the different projections orientations, image them and apply
%colour and transparency
for i=1:num_proj
patch(planes(:,1,i),planes(:,2,i),planes(:,3,i),color(i,:)),alpha(0.9)
end

%Format axis
axis equal

%Plot 'Outgoing X-ray beam'
quiver3(zeros(num_proj,1),zeros(num_proj,1),...
    zeros(num_proj,1),distance*proj_geom.Vectors(:,1),...
    distance*proj_geom.Vectors(:,2),distance*proj_geom.Vectors(:,3));
view(3)
lims=axis;

%Plot 3D axis 
%Define x-axis
plot3(lims(1:2),[0 0],[0 0],'-k') 
%Define y-axis
plot3([0 0],lims(3:4),[0 0],'-k') 
%Define z-axis
plot3([0 0],[0,0],lims(5:6),'-k')
hold off
set(gca,'ZColor','white')
set(gca,'XColor','white')
set(gca,'YColor','white')
end

