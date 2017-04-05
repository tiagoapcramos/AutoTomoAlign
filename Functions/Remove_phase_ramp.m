function [ proj_corr ] = Remove_phase_ramp( this_object,mask )
% REMOVE_PHASE_RAMP Removes phase-ramp from the complex projections
% [ proj_corr ] = Remove_phase_ramp( this_object,mask )
%
% this_object   -   complex 2D projection image
% mask          -   array with same dimensions as 'this_object'. mask=1 for
%                   the pixels corresponding to air in the projection image
%                   and 0 otherwise.
% proj_corr     -   Complex 2D projection image, with respective phase
%                   corrected for possible phase-ramps
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

% Convert mask from logical to numerical array 
mask=double(mask);

% Exclude points not associated with air by converting them to NaN
mask(mask==0)=NaN;

% Isolate phase-ramp (air) from the projection image
this_ramp=angle(this_object).*mask;

% Compute horizontal and vertical gradient field of the phase projection in
% the air sections
[Gx,Gy]=imgradientxy(this_ramp,'centraldifference');

% Creates a rectangular grid with horizontal and vertical pixel coordinates
[X,Y]=meshgrid(1:size(this_ramp,2),1:size(this_ramp,1));

% Remove first order term from the projections
ramp1=median(Gx(~isnan(Gx))).*X+median(Gy(~isnan(Gy))).*Y;
proj_corr=this_object./exp(1i.*ramp1);

% Remove constant term from the projections.
air_level=angle(proj_corr).*mask;
air_level=median(air_level(~isnan(air_level)));
air_level=air_level*ones(size(proj_corr));
proj_corr=proj_corr./exp(1i.*air_level);
end


