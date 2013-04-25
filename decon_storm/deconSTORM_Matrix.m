function [sample_est_mean, sample_est_frames, sample_est_hist, saved_iterations] = ...
   deconSTORM_Matrix(mov, APSF, r, niter, iter_step, alpha, beta, gfactor, fileout, verbose)
% [sample_est_mean, sample_est_frames, sample_est_hist, saved_iterations] = deconSTORM_Matrix(mov, APSF, r, niter, iter_step, alpha, beta, gfactor, fileout, verbose)
%
% In 'matrix' mode, APSF is a rectangular transfer matrix that
%    performs the convolution operation by matrix multiplication.
%
% Inputs:
%  mov - fluorescence movie data; 3d array of size pixels x pixels x frames
%  APSF - point spread function representation (see above).
%  r - background fluorescence intensity
%  niter - number of iterations
%  iter_step - interval between steps which will be output
%  alpha - probability that an active emitter will remain active in the next
%           frame
%  beta - probability that an inactive emitter will become active in the
%           next frame
%  gfactor - gain factor
%  fileout - name of file in which to store output after every batch of iter_step iterations (optional)
%  verbose - 0/1 - whether to print progress reports to the terminal
%
% Outputs:
%  sample_est_mean - Super-resolution estimate of the sample
%  sample_est_frames - Super-resolution estimate of each movie frame
%  sample_est_hist - history of the estimate during the deconvolution process.
%  saved_iterations - vector of the iteration number for each image saved in sample_est_hist
%

% Copyright, 2012
% Eran Mukamel, Hazen Babcock and Xiaowei Zhuang
% contact: eran@post.harvard.edu
%

if nargin<10
   verbose = 1;
end
if nargin<9
   fileout = [];
end

if verbose
   fprintf('********** deconSTORM - Matrix Method ***********\n')
end

[Lx,Ly,nframes] = size(mov);

lx2 = size(APSF,1);
lx = sqrt(lx2); ly = sqrt(lx2);

% normalization constant for each super-resolution pixel
a = sum(APSF,2);

% define temporal filter
filtord = min(floor(log(.0001)/log(alpha)),100);
expfilt2 = alpha.^abs([-filtord:filtord]);
expfilt2 = expfilt2/sum(expfilt2);

saved_iterations = unique(find((mod([1:niter],iter_step)==0) | ([1:niter]==niter)));

a = a * ones(1,nframes);
mov = reshape(mov, Lx*Ly, []);

sample_est = ones(size(a));
sample_est_hist = zeros(lx*ly, length(saved_iterations));

tic
ll=0;
for jiter = 1:niter
   %% deconSTORM iteration
   
   % Estimate of the observed data
   sample_est_down = APSF' * sample_est + r;
   
   EDown = (mov./sample_est_down);
   Econv = APSF * EDown;
   
   sample_est_curr = padarray(sample_est,[0,filtord],0,'pre');
   gam = filter(expfilt2, 1, sample_est_curr, [], 2);
   gam = gam(:,filtord+1:end);
   samplecurr = mean(sample_est,2);
   samplecurr = samplecurr*ones(1,nframes);
   gam = min(gam + beta*samplecurr,1);
   gam = -log(gam) / gfactor;
   
   sample_est = sample_est .* Econv ./ (a+gam);
   
   %% Save intermediate iterations
   if (mod(jiter,iter_step)==0) || (jiter==niter)
      ll=ll+1;
      sample_est_hist(:,ll) = squeeze(mean(sample_est,2));
      tottime = toc;
      if verbose
         fprintf('deconSTORM: iteration %d/%d; time=%3.5g ms/iteration/frame; remaining time=%3.3g s.\n',...
            jiter,niter, 1000*tottime/jiter/nframes, tottime*(niter-jiter)/jiter)
      end
      
      %% Save output
      if ~isempty(fileout)
         sample_est_hist = reshape(sample_est_hist, lx, ly, ceil(niter/iter_step));
         sample_est_frames = reshape(sample_est, lx,ly,nframes);
         sample_est_mean = reshape(mean(sample_est,2),lx,ly);
         
         save(fileout,'sample_est_mean','sample_est_frames','sample_est_hist','saved_iterations');
         if verbose
            fprintf('  Saved current output to %s.\n',fileout)
         end
         
         sample_est_hist = reshape(sample_est_hist, lx*ly, ceil(niter/iter_step));
      end
   end
end

sample_est_hist = reshape(sample_est_hist, lx, ly, ceil(niter/iter_step));
sample_est_frames = reshape(sample_est, lx,ly,nframes);
