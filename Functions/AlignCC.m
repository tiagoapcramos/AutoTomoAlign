function [dis,dis_accum ] = AlignCC( stack_phase_gradient,up_sampling )
% ALIGNCC Align projections using Cross-Correlation
% [dis,dis_accum ] = AlignCC( stack_phase_gradient,up_sampling )
% Uses Cross-correlation measurements to compute translational alignment
% parameters between consecutive stack_phase_gradient.
%
% stack_phase_gradient  -   Gradient intensity field (from the projections)
%                           sorted (by projection angles). Should have 
%                           dimensions [N x M x num_proj]. (The first and 
%                           last projection should be symmetric (0 and 180
%                           degrees for example)
% up_sampling           -   Up-sampling factor during Cross-correlation
%                           measurement. More info in dftregistration
% dis                   -   Matrix with relative displacement vectors
%                           between adjacent projections [u,v]
% dis_accum             -   Matrix with relative displacement vectors in a
%                           global coordinate system (detector's)
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

% Define default up-sampling factor as 1
if nargin<2||isempty(up_sampling),up_sampling=1;end

dis=zeros(size(stack_phase_gradient,3),2);
dis_accum=dis;
h=waitbar(0,'Aligning Projections...');
tic

% Compute displacement vector between consecutive projections
for k=2:size(stack_phase_gradient,3)
d=dftregistration(fft2(stack_phase_gradient(:,:,k)),...
    fft2(stack_phase_gradient(:,:,k-1)),up_sampling);
dis(k,:)=-[d(4) d(3)];
dis_accum(k,:)=sum(dis(1:k,:));
    
time=toc/k;
waitbar(k/size(stack_phase_gradient,3),h,...
    ['Aligning Projections.... Remaining time ',...
    num2str(round((size(stack_phase_gradient,3)-k)*time)),' seconds']);
end
delete(h)

% Computes the centre of rotation (COR) by measuring displacement between
% first and last (flipped) projection. The vertical translation is assumed
% to be due to the accumulation of errors/uncertanties in the measurements
i1=stack_phase_gradient(:,:,1);
iend=fliplr(imtranslate(stack_phase_gradient(:,:,end),dis_accum(end,:)));
overall_d=dftregistration(fft2(i1),fft2(iend),up_sampling);
overall_d=-[overall_d(4) overall_d(3)];

% Distribute vertical uncertanties through all projections and define the
% flobal projections' translations as the cummulative sum of the measured
% displacements between adjacent projections
translatey=linspace(-overall_d(2)/2,overall_d(2)/2,...
    size(stack_phase_gradient,3));
translatex=translatey*0+overall_d(1)/2;
dis_accum=[translatex;-translatey]'+dis_accum;

% Optional - subtracts mean value of vertical alignment parameters
dis_accum(:,2)=dis_accum(:,2)-mean(dis_accum(:,2));
end



