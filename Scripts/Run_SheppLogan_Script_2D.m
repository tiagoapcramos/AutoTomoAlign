%% MAIN DEMO SCRIPT
% Generates a 2D parallel beam phase-contrast tomography problem and solves
% a 2D phase-SIRT reconstruction. The default initialization variables in
% this script return a highly wrapped sinogram which demands the use of a
% non-negativity constrain during the phase-SIRT reconstruction algorithm.
%
% IN THIS DEMO SCRIPT WE USE THE SAME FUNCTIONS AS THOSE FOR A 3D
% SIMULATION SETUP. 
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
% Horizontal dimensions of the 3D phantom
N=200;

% The constant 'max_phase' limits the phantom voxels intensity in order to
% guarantee the maximum admissible phase wrapping in the generated
% projetions.
max_phase=1/10;
% Number of generated projections
num_proj=360;

% GPU properties
% Define maximum available GPU memorry (in Gb) to avoid overload of GPU
% that may lead MATLAB to crash. If you still observe this problem please
% reduce the max_memory value
max_memory=1;

% Number of iterations of the phase-SIRT algorithm for remaining
% tomographic reconstructions
recon_iter=300;
% Define show=1 to visualize the central slice of the reconsructed tomogram
% at a current iteration during reconstruction
recon_options.show=1;
recon_options.max_memory=max_memory;
% Enforce non-negativity constrain during Tomographic Reconstruction
recon_options.non_neg=1;


%% CREATE VOLUME, PROJECTIONS AND PROJECTION GEOMETRY
% Create struct variable 'data' with the defined experimental setup
% Create a 3D Shepp Logan phantom 'volume'
[ data,volume] = CreateTestProblem( N,1,num_proj);
% Pads 'volume' with zeros to accept large translation uncertanties
volume=volume*max_phase;
% volume=volume(:,:,1);
% Updates 'volume' dimensions in 'data'
data.N=size(volume,1);
data.M=1;
% Define matrix with projection parameters
parameters_true=[data.angles_true', data.u_errors, data.v_errors, ...
    data.alpha_errors, data.beta_errors];
%Create
[projection_true,proj_geom]=CreateProjection(volume,parameters_true,...
    max_memory);

%% VISUALIZE PROJECTION GEOMETRY
%!This operation may take some seconds for a large number of projections%
VisualizeProjectionParameters(proj_geom);

%% WRAP SINOGRAM
% Re-order projection_true to create stack of Projections. 
sinogram=permute(projection_true(:,:,1),[3,1,2]);
% WRAP sinogram
sinogram_wrapped=angle(exp(1i.*sinogram));
figure('color','w')
imagesc(squeeze(sinogram_wrapped))
axis equal tight off xy
title('Press any key to continue')
pause

%% TOMOGRAPHIC RECONSTRUCTION WITH phase-SIRT ALGORITHM
X0= Recon_SIRT_wrap_3D(sinogram_wrapped,parameters_true,...
    recon_iter,recon_options);
