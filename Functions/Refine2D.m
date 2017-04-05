function [ parameters_refine ] = Refine2D( volume,stack_projections,...
    parameters,up_sampling,max_memory )
% REFINE2D Refines translational projection parameters
%[ parameters_refine ] = Refine2D( volume,stack_projections,...
%     parameters,up_sampling,max_memory )
% Creates Projections from the current reconstruction and evaluates
% displacement between measured projections and synthetic projections with
% cross-correlation. The cross-correlation alignment is computed with the
% images gradient field intensity
%
% volume              - 3D volume / Tomographic Reconstruction
% stack_projections   - Stack of phase-contrast projections
% parameters          - Matrix with projection parameters for all
%                       projections, in the form [Theta, u, v, alpha, beta]
% up_sampling         - Up-scaling factor for cross-correlation
% max_memory          - Maximum available memory in GPU
% parameters_refine   - Projection parameters after linear parameters
%                       refinement
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

h=waitbar(0,'Creating Synthetic Projections...');

% Create synthetic projections from volume
synth_projections=CreateProjection(volume,parameters,max_memory);
synth_projections=permute(synth_projections,[3 1 2]);

% Compute projections gradient fields
waitbar(0,h,'Creating Gradients');
stack_projections_gradient=phase_gradient(stack_projections);
synth_projections_gradient=phase_gradient(synth_projections);

num_proj=size(stack_projections,3);
parameters_refine=parameters;

%Loop over projections
tic
for k=1:num_proj
   p=stack_projections_gradient(:,:,k);
   p_synth=synth_projections_gradient(:,:,k);
   
   % Calculate displacement
   d=dftregistration(fft2(p),fft2(p_synth),up_sampling);
   d=-[d(4),d(3)];
   
   % Update/Refine projection parameters
   parameters_refine(k,2:3)=parameters_refine(k,2:3)+d;
   time=toc/k;
   waitbar(k/num_proj,h,['Refining Translation Parameters.', ...
       'Remaining time: ',num2str(round((num_proj-k)*time)),' seconds']);
end
delete(h);
end

