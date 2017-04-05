function [ stack_complex,ramps ] = AddPhaseRamps(stack,std0,std1 )
%ADDPHASERAMPS Creates randomly distributed linear backgrounds 
%[ stack_complex,ramps ] = AddPhaseRamps(stack,std0,std1 )
% stack         - Stack of projection images MxNxnum_proj
% stack_complex - Complex stack of projection images with constant
%                 amplitude and phase equal to the input 'stack';
% ramps         - Stack of phase ramps added to the input phase data
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

[M,N,num_proj]=size(stack);

% Create arrays of horizontal and vertical coordinates of the image and
% normalize them from 0 to 1
[X,Y]=meshgrid(1:N,1:M);
X=X/N;
Y=Y/M;
ramps=zeros(size(stack));
rng(3)

% Define random linear (Kxy) and constant (Kz) terms
Kxy=randn(num_proj,2)*std1;
Kz=randn(num_proj,1)*std0/N;

% Creates ramps using plane equation
for k=1:num_proj
   ramps(:,:,k)= Kxy(k,1).*X+Kxy(k,2).*Y+Kz(k);
end

% Create complex stack of projections and add phase ramps
stack_complex=exp(1i.*stack).*exp(1i.*ramps); 

end

