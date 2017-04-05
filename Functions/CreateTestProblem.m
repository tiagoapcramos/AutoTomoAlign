function [ data,volume, vol_geom] = CreateTestProblem( N,M,num_proj,...
    uncert_options)
%CREATETESTPROBLEM Define general simulation setup and phantom
%[ data,volume, vol_geom] = CreateTestProblem( N,M,num_proj,...
%     uncert_options)
% Creates a test problem based on a 3D phantom of size NxNxN
% The number of slices of the phantom (M) can ve controlled by cropping or
% padding the 3D phantom
%
%   N,M   -   Horizontal and Vertical dimensions of the phantom
%   num_proj -  Number of projection angles
%   uncert_options - Struct variable to define errors in the projection
%                   parameters
%        - theta_u - amplitude of random distributed errors in theta
%        - u_u - amplitude of random distributed errors in u
%        - v_u - amplitude of random distributed errors in v
%        - alpha_u - amplitude of random distributed errors in alpha
%        - beta_u - amplitude of random distributed errors in beta
%        - rngs   - integer vector of length 5. Define 'seeds' for
%                   random number generation. Useful for data repeatability
%        - impulse - Bolean variable to include impulse noise in the
%                       translation parameters
%        - noise_perc - density of impulse noise [0 1]
%        - noise_strength - Amplitude of impulse noise
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017


if nargin<4||isempty(uncert_options),uncert_options=struct;end
if nargin<3||isempty(num_proj), num_proj=180;end
if nargin<2||isempty(M),M=N;end

data.N=N;
data.M=M;

% Creates a 3D phantom of dimensions NxNxN
volume=phantom3d_better(data.N);
if M <=N
cut=floor((N-M)/2);
volume=volume(:,:,1+cut:end-cut);
else
pad=round((M-N)/2);
volume=padarray(volume,[0 0 pad]);
end
data.M=size(volume,3);

% Check for default values and missing fields in uncert_options
uncert_options=checkdefaultvalues(uncert_options);

data.num_proj=num_proj;
data.angles_measured=linspace(0,180-180/num_proj,num_proj);

%Define theta errors
rng(uncert_options.rngs(1));
theta_errors = 2*uncert_options.theta_u*(rand(data.num_proj,1) - 0.5);
data.angles_true = data.angles_measured + theta_errors';

%Define u errors
rng(uncert_options.rngs(2));
data.u_errors=2*uncert_options.u_u*(rand(data.num_proj,1)-0.5);
if uncert_options.impulse==true
    rng(uncert_options.rngs(1));
    y_n1=imnoise(zeros(data.num_proj,1),'salt & pepper',...
        uncert_options.noise_perc/100);
    rng(uncert_options.rngs(2));
    y_n2=-imnoise(zeros(data.num_proj,1),'salt & pepper',...
        uncert_options.noise_perc/100);
    data.u_errors=data.u_errors+uncert_options.u_u*(y_n1+y_n2)*...
        uncert_options.noise_strength;
end

%Define v errors
rng(uncert_options.rngs(3));
data.v_errors=2*uncert_options.v_u*(rand(data.num_proj,1)-0.5);
if uncert_options.impulse==true
    rng(uncert_options.rngs(3));
    y_n1=imnoise(zeros(data.num_proj,1),'salt & pepper',...
        uncert_options.noise_perc/100);
    rng(uncert_options.rngs(4));
    y_n2=-imnoise(zeros(data.num_proj,1),'salt & pepper',...
        uncert_options.noise_perc/100);
    data.v_errors=data.v_errors+uncert_options.v_u*(y_n1+y_n2)*...
        uncert_options.noise_strength;
end

%Define alpha errors
rng(uncert_options.rngs(4));                                                                      
data.alpha_errors=2*uncert_options.alpha_u*(rand(data.num_proj,1)-0.5);

%Define beta errors                                                                       
rng(uncert_options.rngs(5));
data.beta_errors=2*uncert_options.beta_u*(rand(data.num_proj,1)-0.5);

if nargout>2
vol_geom=astra_create_vol_geom(size(volume));
end
end

function [uncert_options] = checkdefaultvalues(uncert_options)
% Generates struct variable with default values for uncertanties options in
% case some or all are missing
if ~isfield(uncert_options,'theta_u'), uncert_options.theta_u=0;end
if ~isfield(uncert_options,'u_u'), uncert_options.u_u=0;end
if ~isfield(uncert_options,'v_u'), uncert_options.v_u=0;end
if ~isfield(uncert_options,'alpha_u'), uncert_options.alpha_u=0;end
if ~isfield(uncert_options,'beta_u'), uncert_options.beta_u=0;end
if ~isfield(uncert_options,'rngs'), uncert_options.rngs=[6 3 4 1 2];end
if ~isfield(uncert_options,'impulse'), uncert_options.impulse=false;end
if ~isfield(uncert_options,'noise_perc'), uncert_options.noise_perc=0.1;end
if ~isfield(uncert_options,'noise_strength')
    uncert_options.noise_strength=.5;
end
end
