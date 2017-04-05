function [ Gmag,Gy ] = phase_gradient( P,fd_method )
% PHASE_GRADIENT Computes phase projections gradient field
% [ Gmag,Gy ] = phase_gradient( P,fd_method )
% Computes the image gradient field, assuming that the gradients in x and y
% direction are within ]-pi pi[
% If the number of arguments in the output is equal to 1, the function
% returns the magnitude of the phase gradient field. If the number of
% output arguments is equal to 2, the function returns both the gradient in
% x (Gx) and y (Gy)
%
% P         - Stack of Phase projections - dimension [M,N,num_proj]
% fd_method - Method for finite differences computation = 'forward' or
%             'central'
% Gmag      - Magnitude of the phase gradient field (1 output variable) or
%             Phase gradient field in the x direction (2 output variables)
% Gy        - Phase gradient field in the y direction
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

% Set default finite difference method to 'forward'
if nargin<3||isempty(fd_method),fd_method='forward';end

% Computes Gradient field in x and y direction through a finite differences
% scheme. The operation is made in the complex domain to allow for phase
% wrapping in the projections
if strcmp(fd_method,'forward')==1
    Gx=angle(exp(1i.*P(:,2:end,:))./exp(1i.*P(:,1:end-1,:)));
    Gx(:,end+1,:)=Gx(:,end,:);
    Gy=angle(exp(1i.*P(2:end,:,:))./exp(1i.*P(1:end-1,:,:)));
    Gy(end+1,:,:)=Gy(end,:,:);
elseif strcmp(fd_method,'central')==1
    Gx=angle(exp(1i.*P(:,3:end,:))./exp(1i.*P(:,1:end-2,:)))/2;
    Gx=padarray(Gx,[0 1 0],'replicate','both');
    Gy=angle(exp(1i.*P(3:end,:,:))./exp(1i.*P(1:end-2,:,:)))/2;
    Gy=padarray(Gy,[1 0 0],'replicate','both');
end

if nargout>1 , Gmag=Gx;
else Gmag=sqrt(Gx.^2+Gy.^2);
end
end

