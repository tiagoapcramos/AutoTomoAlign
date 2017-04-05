function [ X,parameters] = Refine3D( stack_sinogram,parameters0,...
    X,max_out_iter,jac_options,recon_iter,recon_options,...
    lsq_options,max_memory )
% REFINE3D Refines projection parameters using an iterative solver
%[ X,parameters] = Refine3D( stack_sinogram,parameters0,X,max_out_iter,...
%     jac_options,recon_iter,recon_options,lsq_options,max_memory )
% Estimates a new set of projection parameters by solving a non-linear
% least squares problem.
%     
% stack_sinogram   -   stack of sinograms. have dimensions [N,M,theta]
% parameters0      -   matrix with projection parameters
% X (input)        -   Current reconstructed volume
% X (output)       -   Reconstructed volume after parameters optimisation
% max_out_iter     -   Number of outer iterations in the optimisation
% jac_options      -   Jacobian options
%                      jac_options.step_size - [1x5] vector with jacobian
%                                              stepsizes
%                      jac_options.method - 'forward' or 'central' for 
%                                           finite differences scheme
% recon_iter      -   number of iterations to run
% recon_options   -   .non_neg - non negativity constrain. may help in the
%                        presence of large wrapping areas
%                        .show=1 shows image evolution
%                        .slice  slice to monitor while reconstructing
%                        .max_memory - maximum memory of GPU card available
% lsq_options     - optimisation options check 'doc optimoptions' or
%                   'doc lsqnonlin'
% parameters      -   matrix with projection parameters after optimisation
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

parameters=parameters0;
for k=1:max_out_iter
    
    % Define Cost Function
    f_fun=@(p) Cost_Function(permute(stack_sinogram,[2,3,1]),X,p,...
        jac_options.method,jac_options.step_size,max_memory );
    
    % Run Levenberg-Marquardt Optimisation
    parameters=lsqnonlin(f_fun,parameters,[],[],lsq_options);
    
    % Make new tomographic reconstruction
    X= Recon_SIRT_wrap_3D(stack_sinogram,parameters,recon_iter,...
        recon_options );
end
end
    