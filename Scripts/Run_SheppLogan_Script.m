%% MAIN DEMO SCRIPT
% Generates a 3D parallel beam phase-contrast tomography problem with
% translational and angular uncertanties in the projections parameters.
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

%% INITIALIZATION VARIABLES - Please define
% Include the paths for the Functions and External folders from the same
% release as this script file
functions_path=fullfile('..','Functions');
external_path=fullfile('..','External');
% Include the paths for the mex\ and tools\ folders from your ASTRA toolbox
ASTRA_path={'C:\ASTRA\astra-1.6\mex','C:\ASTRA\astra-1.6\tools'};
addpath(functions_path)
addpath(external_path)
addpath(ASTRA_path{:})

% Volume and projection parameters generation
% Horizontal dimensions of the 3D phantom (before padding)
N=300;
% Vertical dimension of the 3D phantom (before padding)
M=300;

% When simulating large translational uncertanties, in order to guarantee
% that all projections are inside the field of view, pad the generated 
% phantom with zeros.
padxy=20;
padz=20;

% The constant 'max_phase' limits the phantom voxels intensity in order to
% guarantee the maximum admissible phase wrapping in the generated
% projetions.
max_phase=1/20;
% Number of generated projections
num_proj=360;

% Uncertanties/Errors amplitude (pixels for u and v, and degrees for theta,
% alpha and beta)in the projection parameters
uncert_options.theta_u=1;
uncert_options.u_u=10;
uncert_options.v_u=10;
uncert_options.alpha_u=5;
uncert_options.beta_u=5;
uncert_options.rngs=[6 3 4 1 2];
% To include additional 'impulse-like' noise in the translational
% projection parameters set uncert_options.impulse to true
uncert_options.impulse=false;
% noise_perc controls the frequency of the impulse noise uncertanties
uncert_options.noise_perc=0.1;
% noise_strength defines the amplitude of the noise uncertanties
uncert_options.noise_strength=.5;
% Standard deviation for the 1st order term (linear) in phase-ramps
uncert_options.ramps_std1=0.5;
% Standard deviation for the 0th order term (constant) in phase-ramps
uncert_options.ramps_std0=0.3;

% GPU properties
% Define maximum available GPU memorry (in Gb) to avoid overload of GPU
% that may lead MATLAB to crash. If you still observe this problem please
% reduce the max_memory value
max_memory=1;

% Gradient Options - Finite differences
fd_method='forward';

% Cross-Correlation alignment up sampling factor
up_sampling=100;

% Reconstruction Options
% Number of iterations of the phase-SIRT algorithm for the auxiliar
% reconstruction (during Phase Ramp Removal)
recon_iter_aux=10;
% Number of iterations of the phase-SIRT algorithm for remaining
% tomographic reconstructions
recon_iter=300;
% Define show=1 to visualize the central slice of the reconsructed tomogram
% at a current iteration during reconstruction
recon_options.show=1;
recon_options.max_memory=max_memory;
% Enforce non-negativity constrain during Tomographic Reconstruction
recon_options.non_neg=false;

% Segmentation Options (during Phase Ramp Removal)
% chose 'slice' for a slice-by-slice binary segmentation. Chose 'volume'
% instead for a 3D segmentation.
seg_method='volume';

% Optimisation options
% Maximum number of 'outer' iterations
max_out_iter=5;
% Optimization options - please check 'doc optimset'
tol_fun = 1e-8; % optimize - tolerance [default 1E-8]
tol_x = 1e-8; % optimize - tolerance [default 1E-8]
max_funevals = 2000; %optimize - max iter [default 1E4]
max_in_iter = 15; %optimize - max iter [default 1E4]
lambda0=1; %Initial Lagrange multiplier
lsq_options = optimset('Algorithm',{'levenberg-marquardt',lambda0},...
   'TolFun',tol_fun,'TolX',tol_x,'MaxFunEvals',max_funevals,...
   'MaxIter',max_in_iter','Display','iter-detailed','Jacobian','on');
% Jacobian finite differences method. Chose between 'forward' and 'central'
jac_options.method='forward';
% Finite differences step-size for Jacobian calculation
jac_options.step_size=[1,1,1,1,1]*1e-3   ;
clear tol_fun tol_x max_funevals max_iter

%% CREATE VOLUME, PROJECTIONS AND PROJECTION GEOMETRY
% Create struct variable 'data' with the defined experimental setup
% Create a 3D Shepp Logan phantom 'volume'
[ data,volume] = CreateTestProblem( N,M,num_proj,uncert_options);

% Pads 'volume' with zeros to accept large translation uncertanties
volume=padarray(volume,[padxy padxy padz])*max_phase;

% Updates 'volume' dimensions in 'data'
data.N=size(volume,1);
data.M=size(volume,3);

% Define matrix with projection parameters
parameters_true=[data.angles_true', data.u_errors, data.v_errors, ...
    data.alpha_errors, data.beta_errors];

%Create projections
[projection_true,proj_geom]=CreateProjection(volume,parameters_true,...
    max_memory);

%% VISUALIZE PROJECTION GEOMETRY
%!This operation may take some seconds for a large number of projections!%
VisualizeProjectionParameters(proj_geom);

%% ADD PHASE-RAMPS AND CREATE STACK OF COMPLEX PROJECTIONS
% Re-order projection_true to create stack of Projections. 
stack_object=permute(projection_true,[3,1,2]);

% Add phase ramps and convert to complex projection images of constant
% amplitude.
[ stack_object,ramps0 ] = AddPhaseRamps(stack_object, ...
    uncert_options.ramps_std0,uncert_options.ramps_std1 );

%% VISUALIZE PROJECTIONS
ViewProjections(stack_object);
%% TOMOGRAPHIC RECONSTRUCTION BEFORE ALIGNMENT
% Define matrix with measured projection parameters
parameters0=zeros(num_proj,5);
parameters0(:,1)=data.angles_measured;
X0= Recon_SIRT_wrap_3D(angle(stack_object),parameters0,...
    recon_iter,recon_options);


%% CHECK HISTOGRAM OF DIFFERENCE BETWEEN REAL PHASE AND WRAPPED PHASE
figure('color','w')
check_lims=angle(stack_object)-unwrap(angle(stack_object));
histogram(check_lims);
if max(check_lims(:))>pi || min(check_lims(:))<-pi
    title('The amount of wrapping or phase ramp may be too high')
else
    title('Required conditions for maximum wrapping confirmed')
end
%% ALIGN PROJECTIONS USING CROSS-CORRELATION IN THE GRADIENT DOMAIN
% Compute gradient field intensity of the phase projection data
stack_phase_gradient=phase_gradient(angle(stack_object),fd_method);
% Performs Cross-Correlation alignment
[~,dis_accum ] = AlignCC( sqrt(stack_phase_gradient),up_sampling);

%% VISUALIZE PROJECTIONS
% Check the results from the CC alignment
ViewProjections(stack_object,dis_accum);
%% CREATE AUXILIAR RECONSTRUCTION FOR PHASE RAMP REMOVAL
% Define matrix with projection parameters including the projections 
% relative translations determined by Cross-Correlation
parameters=zeros(num_proj,5);
parameters(:,1)=data.angles_measured;
parameters(:,2:3)=dis_accum;

% Creates a circular mask
mask=repmat(Mask_Circle(data.N-2,data.N),[1 1 data.M]);

%Auxiliar Tomographic Reconstruction
rec_aux= Recon_SIRT_wrap_3D(angle(stack_object),parameters,...
    recon_iter_aux,recon_options,mask);

%% SEGMENT AUXILIAR RECONSTRUCTION AND CREATE SYNTHETIC PROJECTIONS
% Segment auliar reconstruction
rec_aux_segmented = SegmentVolume(rec_aux,seg_method);

% Creates Synthetic Projections
proj_auxiliar=CreateProjection(rec_aux_segmented,parameters,max_memory);
stack_aux=permute(proj_auxiliar,[3 1 2]);

% Define stack_aux as logical array with 1 for regions correspoing to air
% and 0 otherwise
stack_aux=stack_aux==0;
ViewProjections(double(stack_aux),dis_accum);
%% REMOVE PHASE RAMP
% Corrects for phase-ramps in the complex stack of projections by using 
% the 'air' values, assigned by stack_aux, as reference.
stack_corr=zeros(size(stack_object));
stack_aux=logical(stack_aux);
warning('off','all')
for k=1:num_proj
    this_object=stack_object(:,:,k);
    this_mask=stack_aux(:,:,k);
    object_corr=Remove_phase_ramp(this_object,this_mask);
    stack_corr(:,:,k)=object_corr;
    disp(['Done ',num2str(k),'/',num2str(num_proj)]);
end
warning('on','all')
%% VISUALIZE PROJECTIONS
% Check the results from Phase Ramp Removal
ViewProjections(stack_corr,dis_accum);

%% ALIGN CROSS CORRELATION (again - optional!)
stack_phase_gradient=sqrt(phase_gradient(angle(stack_corr),fd_method));
[dis,dis_accum ] = AlignCC( stack_phase_gradient,up_sampling);
parameters(:,2:3)=dis_accum;

%% TOMOGRAPHIC RECONSTRUCTION 
% Creates a tomographic reconstruction after cross-correlation alignment 
% and Phase Ramp Removal
X_CTA= Recon_SIRT_wrap_3D(angle(stack_corr),parameters,recon_iter,...
    recon_options,mask);
%% MAKE 2D REFINEMENT WITH Cross-Correlation on synthetic projections
% (Optional)
% parameters_refine=Refine2D(X_CTA,angle(stack_corr),parameters,...
%     up_sampling,max_memory);
%% MAKE 3D REFINEMENT WITH LEVENBERG-MARQUARDT OPTIMISATION
% ! This operation can be time-consuming !
parameters_before_FPP=parameters;
[ X_FPP,parameters_FPP] = Refine3D(angle(stack_corr),...
    parameters_before_FPP,X_CTA,max_out_iter,jac_options,recon_iter,...
    recon_options,lsq_options,max_memory );
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%       END OF ALIGNMENT AND TOMOGRAPHIC RECONSTRUCTION       %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% COMPUTES COST-FUNCTION
% Computes Cost-Function using the projection parameters at different
% stages of the proposed alignment algorithm
p0=permute(CreateProjection(X0,parameters0,max_memory),[3,1,2]);
r0=angle(exp(1i.*angle(stack_object)-1i.*p0));
r0=0.5*sqrt(sum(r0(:).^2));

X1=Recon_SIRT_wrap_3D(angle(stack_object),parameters,recon_iter,...
    recon_options);
p1=permute(CreateProjection(X1,parameters,max_memory),[3,1,2]);
r1=angle(exp(1i.*angle(stack_object)-1i.*p1));
r1=0.5*sqrt(sum(r1(:).^2));

p2=permute(CreateProjection(X_CTA,parameters,max_memory),[3,1,2]);
r2=angle(exp(1i.*angle(stack_corr)-1i.*p2));
r2=0.5*sqrt(sum(r2(:).^2));

p3=permute(CreateProjection(X_FPP,parameters_FPP,max_memory),[3,1,2]);
r3=angle(exp(1i.*angle(stack_corr)-1i.*p3));
r3=0.5*sqrt(sum(r3(:).^2));

%% SAVES RESULTS
clear projection_true rec_aux rec_aux_segmented this_mask this_object ...
      p0 p1 p2 p3
save(['results_',date],'-v7.3');
