function [ r0,J ] = Cost_Function(projection_true,X,parameters, method,...
                    step_size,max_memory )
%COST_FUNCTION Defines Cost-funciton for optimisation problem
% [ r0,J ] = Cost_Function(projection_true,X,parameters, method,...
%                     step_size,max_memory )
% Computes residual vector and Jacobian of the optimization problem. 
% For more info see 'doc lsqnonlin'
%
% projection_true - stack of phase projections ASTRA convention [N,theta,M]
% X               - 3D volume/ Tomographic Reconstruction
% parameters      - Matrix with projection parameters for all
%                   projections, in the form [Theta, u, v, alpha, beta]
% method          - Jacobian finite differences method. Chose between 
%                   'forward' and 'central'
% step_size       - Finite differences step-size for Jacobian calculation
% max_memory      - Maximum available memory in GPU
% r0              - Residual Vector
% J               - Jacobian Matrix
%
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

% Define default step-size and finite differences method for the Jacobian
% matrix calculation
if nargin<5||isempty(step_size),step_size=[0.1,0.1,0.1,0.1,0.1]*1e-3;end
if nargin<4||isempty(method),method='forward';end

num_proj=size(parameters,1);
proj_size=numel(projection_true)/num_proj;
m=numel(projection_true);
n=numel(parameters);
y0=permute(CreateProjection(X,parameters,max_memory),[3,1,2]);

% COMPUTE JACOBIAN MATRIX IF REQUIRED (2 arguments in the function output)
% Allocate required memory for Jacobian matrix
if nargout>1
    J=spalloc(m,n,5*numel(projection_true));

% Compute Jacobian using forward differences scheme
if strcmp(method,'forward')==1
   % Loop over different parameters: theta,u,v,alpha,beta
   for l=1:5        
       aux=zeros(1,5);aux(l)=1;
       delta_param=repmat(step_size.*aux,[num_proj,1]);
       
       % Compute projections
       % Permute is necessary to order the residual vector with Tiago's
       % convention. 
       y2=permute(CreateProjection(X,parameters+delta_param,...
           max_memory),[3,1,2]);
       
       %Compute finite difference
       delta_r=y2(:)-y0(:);
       col_indexes=reshape(repmat((1:num_proj)+(l-1)*...
           num_proj,proj_size,1),1,[]);
       J=J+sparse(1:numel(y2),col_indexes,delta_r./step_size(l),m,n);
       display('Jacobian computed')
   end
   
% Compute Jacobian using central differences scheme
elseif strcmp(method,'central')==1
   %Loop over different parameters: Theta, t_horiz, t_vert,alpha,beta
   for l=1:5       
       aux=zeros(1,5);aux(l)=1;
       delta_param=repmat(step_size.*aux,[num_proj,1]);
       
       % Compute projections
       % Permute is necessary to order the residual vector with Tiago's
       % convention. 
       y1=permute(CreateProjection(X,parameters-delta_param,...
           max_memory),[3,1,2]); %y-1
       y2=permute(CreateProjection(X,parameters+delta_param,...
           max_memory),[3,1,2]); %y+1

       %Compute finite difference
       delta_r=(y2(:)-y1(:))/2;
       col_indexes=reshape(repmat((1:num_proj)+(l-1)*num_proj,...
           proj_size,1),1,[]);
       J=J+sparse(1:numel(y2),col_indexes,delta_r./step_size(l),m,n);
       display('Jacobian computed')
        
   end
 
end
end

% Compute Residual Vector
r0=angle(exp(1i.*y0)./exp(1i.*permute(projection_true,[3,1,2])));
r0=r0(:);
display('residuals computed')
end