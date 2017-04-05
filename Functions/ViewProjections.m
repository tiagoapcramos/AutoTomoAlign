function [ h ] = ViewProjections( stack,trans )
% VIEWPROJECTIONS displays animation of the stack of projections
%[ h ] = ViewProjections( stack,trans )
% Images projections in 'stack' applying  (if existing) the computed
% translations in 'trans'. The projections can be either. complex or
% double. In the case of complex projections, this function images the
% projections phase values.
%
% stack    - stack of projections (sorted by projection angle)
% trans    - Matrix with the linear projection parameters [u,v] for all
% projections
% if 'stack' is complex then the function will image the phase of the array
%
% case_n 0 - stack is complex no translation
% case_n 1 - stack is complex with translations
% case_n 2 - stack is real (phase) no translation
% case_n 3 - stack is real (phase) with translation
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

% Check for type of input variable
if (nargin<2||isempty(trans))&&(isreal(stack)==0),case_n=0;
elseif ~(nargin<2||isempty(trans))&&(isreal(stack)==0),case_n=1;
elseif (nargin<2||isempty(trans))&&(isreal(stack)==1),case_n=2;
elseif ~(nargin<2||isempty(trans))&&(isreal(stack)==1),case_n=3;
else error('Please check input variables')
end
                                                                           
h=figure('color','w');

% Loop over the different projections and image. The default pausing time
% between projections was set to 0.001 seconds
for k=1:size(stack,3)
switch case_n
    case 0,imagesc(angle(stack(:,:,k)))
    case 1,imagesc(angle(imtranslate(stack(:,:,k),trans(k,:))))
    case 2,imagesc(stack(:,:,k))
    case 3,imagesc(imtranslate(stack(:,:,k),trans(k,:)));
end
axis equal tight off xy
title([num2str(k),'/',num2str(size(stack,3))])
pause(0.001)
end
end

