function [ recon,residuals ] = Recon_SIRT_wrap_3D(stack_sinogram,...
    parameters,n_iters,options,mask )
% RECON_SIRT_WRAP_3D Makes 3D parallel beam tomographic reconstruction
%[ recon,residuals ] = Recon_SIRT_wrap_3D(stack_sinogram,...
%     parameters,n_iters,options,mask )
% 3D iterative tomographic reconstruction algorithm for phase-contrast
% projection data. The modifications made to the SIRT algorithm allow the
% presence of wrapping in the phase projections within reasonable
% limitations (for real data applications)
%
% stack_sinogram   -   stack of sinograms. have dimensions [N,M,theta]
% parameters -   matrix with projection parameters : [theta,u,v,alpha,beta]
% n_iters          -   number of iterations to run
% options          -   .non_neg - non negativity constrain. may help in the
                       %presence of large wrapping areas
                       %.show=1 shows image evolution
                       %.slice  slice to monitor while reconstructing
                       %.max_memory - maximum memory of GPU card available
% mask             -   3D 'Support' mask to multiply to Reconstructed
                        % volume if desired
% recon            -   3D Tomographic Reconstruction
% residuals        -   Vector with Cost-function at each iteration
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

% Define default number of iterations and reconstruction options
if nargin<3||isempty(n_iters),n_iters=100;end
if nargin<4 || isempty(options)
    options.non_neg=1;
    options.show=1;
    options.max_memory=10;
end
if nargout>1,residuals=zeros(n_iters,1);end

% Re-order stack of Projections to respect ASTRA's convention
if size(stack_sinogram,3)==1 %input is 2D and a single sinogram
    stack_sinogram=permute(stack_sinogram,[3,1,2]);
end

% Create ASTRA volume geometry
N=size(stack_sinogram,2);
M=size(stack_sinogram,1);
recon=zeros(N,N,M);
vol_geom=astra_create_vol_geom(N,N,M);

%Create ASTRA projection geometry according to the projection parameters
proj_geom=astra_create_proj_geom('parallel3d',1,1,M,...
    N,parameters(:,1)*pi/180);
proj_geom=astra_geom_2vec(proj_geom);
u=zeros(size(parameters,1),3);
v=u;
for i=1:size(parameters,1)
    R=rotz(parameters(i,1))*rotx(parameters(i,4))*roty(parameters(i,5));
    u(i,:)=(R*[1 0 0]')';                      
    v(i,:)=(R*[0 0 1]')';   
end
                    
proj_geom.Vectors(:,7:12)=[u v];
proj_geom.Vectors(:,4:5)=proj_geom.Vectors(:,7:8).*[parameters(:,2)...
    parameters(:,2)];
proj_geom.Vectors(:,6)=proj_geom.Vectors(:,12).*parameters(:,3);
proj_geom.Vectors(:,1:3)=cross(proj_geom.Vectors(:,7:9)',...
    proj_geom.Vectors(:,10:12)')'; 

% SIRT initialization and creation of 'matrices' C and R (large and
% diagonal = unecessary, similar to .* matrix multiplication
projection_true=permute(stack_sinogram,[2,3,1]);

%Define forward and back funcction operators
Afun=@(X) forwardproj(X,proj_geom,vol_geom,options);
Atfun=@(proj) backproj(proj,proj_geom,vol_geom);

C=Atfun(ones(size(projection_true)));
C=1./(C);C(isinf(C))=0;
R=Afun(ones(size(recon)));
R=1./(R);R(isinf(R))=0;
% end initialization %

if options.show==1
    figure('color','w');
end
if (options.show==1||options.show==2)&&(~isfield(options,'slice'))
    options.slice=round(M/2);
end

%%% Start Iteration %%%
aux=exp(1i.*projection_true); %to avoid repeat this operation
for k=1:n_iters

    % Calculate residual vector
    r=angle(aux./exp(1i.*Afun(recon)));

    % Calculate Cost-Function if required as output
    if nargout>1,residuals(k)=0.5*r(:)'*r(:);end

    % Update reconstruction
    recon=recon+C.*(Atfun(R.*r));

    % Apply 'mask support' if required
    if nargin>4,recon=recon.*mask;end

    % Enforce non-negativity if required
    if options.non_neg==true,recon=max(0,recon);end

    % Show 1 slice of the reconstructed volume at each iteration during
    % reconstruction
    if options.show==1
        if nargout>1
            figure(1)
            subplot(121),imagesc(recon(:,:,options.slice)),
            axis equal tight off xy,title(num2str(k))
            subplot(122),plot(residuals(1:k)),xlabel('iteration number'),
            ylabel('sum of sq. residuals')
            drawnow
        else
            figure(1)
            imagesc(recon(:,:,options.slice)),axis equal tight off xy,
            title(num2str(k)),colorbar
           drawnow
        end
    end
end

function [proj] = forwardproj(X,proj_geom,vol_geom,options)
% Forward Projection operator. Generates tomographic projections from the
% volume X, according to the projector geometry defined by proj_geom. This
% operation is equivalent to A*x
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

N=size(X,1);M=size(X,3);
n_proj=size(proj_geom.Vectors,1);
proj=zeros(N,n_proj,M);

% Calculate number of bytes required for the GPU to storage the volume X
% and generated projections
req_memory=(numel(proj)+N*N*M)*8/(1024^3); %In Gbytes (double precision)

if req_memory<options.max_memory
    [~,proj]=astra_create_sino3d_cuda(X,proj_geom,vol_geom);
    astra_mex_data3d('clear');
else
    
% If the available GPU memory is not sufficient to generate all
% projections at once, the function automatically generates smaller
% stacks of proejctions that are merged into a single array later.
div=ceil(req_memory/options.max_memory);
n_partial_proj=ceil(n_proj/div);
indexes=1:n_partial_proj;
for k=1:div
    proj_geom_k=proj_geom;
    proj_geom_k.Vectors=proj_geom.Vectors(indexes,:);
    [~,proj_k]=astra_create_sino3d_cuda(X,proj_geom_k,vol_geom);
    proj(:,indexes,:)=proj_k;
    indexes=indexes+n_partial_proj;
    indexes(indexes>n_proj)=[];
    astra_mex_data3d('clear');
end
end

function [BP] = backproj(proj,proj_geom,vol_geom)
% Back Projection operator. Performs a Back-Projection operation from the
% projections proj, according to the projector geometry defined by 
% proj_geom. This operation is equivalent to A'*x
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

projection_id=astra_mex_data3d('create','-sino',proj_geom,proj);
BP_id=astra_mex_data3d('create','-vol',vol_geom,0);
cfg=astra_struct('BP3D_CUDA');
cfg.ReconstructionDataId=BP_id;
cfg.ProjectionDataId=projection_id;
alg_id=astra_mex_algorithm('create',cfg);
astra_mex_algorithm('iterate',alg_id);
BP=astra_mex_data3d('get',BP_id);
astra_mex_algorithm('clear')
astra_mex_data3d('clear')

