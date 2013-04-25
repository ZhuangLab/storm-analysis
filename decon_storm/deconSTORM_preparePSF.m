function PSF = deconSTORM_preparePSF(sigma,Npixels,npixels)
% PSF = deconSTORM_preparePSF(sigma,Npixels,npixels)
%
% Create a Gaussian point spread function (PSF) matrix
%
% Inputs:
% sigma - Standard deviation parameter of the Gaussian PSF shape, in pixels
% Npixels - linear dimension, in pixels, of the image field of view
% (assumed square)
% npixels - Linear dimension (in pixels) of the super-resolution estimate
%
% Output:
% PSF - point spread function, of size npixels x npixels
% APSF - transfer function mapping from between the high- and
% low-resolution images, of size (npixels^2) x (Npixels^2)
%

% Copyright, 2012
% Eran Mukamel, Hazen Babcock and Xiaowei Zhuang
% Contact: eran@post.harvard.edu
%

dsamp = npixels/Npixels;

if mod(npixels,Npixels)~=0
   error('npixels must be an integer multiple of Npixels.')
end

xvec = ([1:npixels]-1/2)*Npixels/npixels;
Xvec = ([1:Npixels]-1/2)*Npixels/Npixels;

[xx,yy] = meshgrid(xvec);
PSF = exp(-((xx-(Npixels+1)/2).^2 + (yy-(Npixels+1)/2).^2)/(2*sigma^2));

if dsamp>1
   dsampvec = [dsamp/2:dsamp:npixels];
elseif dsamp==1
   dsampvec = [1:npixels];
end

norm = sum(sum(PSF(dsampvec,dsampvec)));
PSF = PSF / norm;

