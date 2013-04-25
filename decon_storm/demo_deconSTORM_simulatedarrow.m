%
% Demo script for deconSTORM, using simulated fluorescence data.
%

% Copyright, 2012
% Eran Mukamel, Hazen Babcock and Xiaowei Zhuang
% contact: eran@post.harvard.edu
%

%% Load sample STORM movie data

% Fluorescence movie data from imaging of microtubules:
load('deconSTORM_simulatedarrowdata.mat','mov')

%% Set parameters

% Npixels is the linear dimension, in pixels, of the image field of view (assumed square)
% nframes is the number of movie frames
[Npixels,~,nframes] = size(mov);

% Factor by which we will sub-sample each image dimension to create a
% super-resolution sample estimate
dsamp = 8; 

% Linear dimension (in pixels) of the super-resolution estimate
npixels = Npixels*dsamp;

% Probability that an active emitter remains active in the next frame
alpha = 1/2;

% Probability that an inactive emitter will become active in the next frame
beta = 1/120;

% Background fluorescence intensity, in photons per pixel per frame
r =1;

% Standard deviation parameter of the Gaussian PSF shape, in pixels
sigma = 1; 

% Gain parameter for deconSTORM
gfactor = 256;

%% Generate the point spread function (PSF)

% PSF is a Gaussian point spread function, of size npixels x npixels.
% APSF is a matrix of size (npixels^2) x (Npixels^2). Each column of APSF
% is the expected image (flattened into a column vector) of an emitter at a
% single point in the super-resolution field of view.


%% ---------- METHOD 1: Run deconSTORM using Matrix method
% This method may be computationally faster when there is sufficient memory
% to compute the transfer matrix, APSF. Also, this method does not assume
% periodic boundary conditions.

% Generate the point spread function transfer matrix, APSF
[APSF] = deconSTORM_prepareAPSF(sigma,Npixels,npixels);

% Number of iterations of deconSTORM
niter = 1000;

% Interval between iterations at which to report the output
iter_step = 50;

% Name of file in which to store results
fileout = 'deconSTORM_microtubule_results.mat';

verbose = 1; 

[sample_est_mean, sample_est_frames, sample_est_hist, saved_iterations] = deconSTORM_Matrix(mov, APSF, ...
   r, niter, iter_step, alpha, beta, gfactor, fileout, verbose);

%% ---------- METHOD 1: Run deconSTORM using Convolution method
% This method may be preferrable for large images, for which the transfer
% matrix APSF is too large to store in memory.  This method assumes
% periodic boundary conditions for the image.

% Generate the point spread function, PSF
[PSF] = deconSTORM_preparePSF(sigma,Npixels,npixels);

% Number of iterations of deconSTORM
niter = 1000;

% Interval between iterations at which to report the output
iter_step = 50;

% Name of file in which to store results
fileout = 'deconSTORM_microtubule_results.mat';
verbose = 1;

[sample_est_mean, sample_est_frames, sample_est_hist] = deconSTORM_Conv(mov, PSF, ...
   r, niter, iter_step, alpha, beta, gfactor, fileout, verbose);

