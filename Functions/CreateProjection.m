function [ projection,proj_geom ] = CreateProjection( volume,...
    parameters,max_memory )
%CREATEPROJECTION Creates tomographic projections
% [ projection,proj_geom ] = CreateProjection( volume,...
%     parameters,max_memory )
% Creates a 3D matrix with a stack of sinograms by computing discrete ...
%integrations over 'volume' using the projection parameters in 'parameters'
%
% volume              - 3D volume / Tomographic Reconstruction
% parameters          - Matrix with projection parameters for all
%                       projections, in the form [Theta, u, v, alpha, beta]
%                       (:,1) - projection angles (degrees)
%                       (:,2) - horizontal granslations (Pixels)
%                       (:,3) - vertical translations (Pixels)
%                       (:,4) - alpha (Degrees)
%                       (:,5) - beta (Degrees)
% max_memory          - Maximum available memory in GPU
% projection          - 3D matrix with stack of sinograms (ASTRA's
%                       convention)
% proj_geom           - Structure array with 3D parallel beam geometry (see
%                       ASTRA's toolbox documentation)
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

if size(parameters,2)<5
    aux=zeros(size(parameters,1),5);
    aux(:,1:size(parameters,2))=parameters;
    parameters=aux;
end

%Create volume Geometry
vol_geom=astra_create_vol_geom(size(volume));

%Create 3D parallel beam geometry
proj_geom=astra_create_proj_geom('parallel3d',1,1,size(volume,3),...
    size(volume,1),parameters(:,1)*pi/180);
proj_geom=astra_geom_2vec(proj_geom);
u=zeros(size(parameters,1),3);
v=u;

% Define 3D Transformation/rotation matrix
for i=1:size(parameters,1)
    R=rotz(parameters(i,1))*rotx(parameters(i,4))*roty(parameters(i,5));
    u(i,:)=(R*[1 0 0]')';                     
    v(i,:)=(R*[0 0 1]')';   
end                    
proj_geom.Vectors(:,7:12)=[u v];   
proj_geom.Vectors(:,4:5)=proj_geom.Vectors(:,7:8).*...
    [parameters(:,2) parameters(:,2)]; 
proj_geom.Vectors(:,6)=proj_geom.Vectors(:,12).*parameters(:,3);
proj_geom.Vectors(:,1:3)=cross(proj_geom.Vectors(:,7:9)',...
    proj_geom.Vectors(:,10:12)')';  

% Check for available memory and dvide projections in subsets if required
[~,N,M]=size(volume);
num_proj=size(parameters,1);
req_memory=M*N*num_proj*8/(1024^3);%In Gbytes (double precision)
if req_memory<max_memory    
    [~,projection]=astra_create_sino3d_cuda(volume,proj_geom,vol_geom);
else
    
% If the available GPU memory is not sufficient to generate all
% projections at once, the function automatically generates smaller
% stacks of proejctions that are merged into a single array later.
projection=zeros(N,num_proj,M);
div=ceil(req_memory/max_memory);
n_partial_proj=ceil(num_proj/div);
indexes=1:n_partial_proj;
for k=1:div
    proj_geom_k=proj_geom;
    proj_geom_k.Vectors=proj_geom_k.Vectors(indexes,:);
    [~,aux]=astra_create_sino3d_cuda(volume,proj_geom_k,vol_geom);
    projection(:,indexes,:)=aux;
    indexes=indexes+n_partial_proj;
    indexes(indexes>num_proj)=[];
    astra_mex_data3d('clear');
end
    astra_mex_data3d('clear');
end