function [sample_est_mean, sample_est_frames, sample_est_hist, saved_iterations] = ...
   deconSTORM_Conv(mov, APSF, r, niter, iter_step, alpha, beta, gfactor, fileout, verbose)
% [sample_est_mean, sample_est_frames, sample_est_hist, saved_iterations] = deconSTORM_Conv(mov, APSF, r, niter, iter_step, alpha, beta, gfactor, fileout, verbose)
%
%    In 'convolution' mode, APSF is a 2d matrix, the same size as the
%    super-resolution image, which will be convolved with the image
%    estimate at each iteration.
%
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
   fprintf('********** deconSTORM - Convolution Method ***********\n')
end

[Lx,Ly,nframes] = size(mov);

[lx,ly] = size(APSF);

PSFk = fft2(APSF);
PSFk = repmat(PSFk,[1,1,nframes]);

dsampx = lx/Lx;
dsampy = ly/Ly;

% normalization constant
a = sum(APSF(:))/(dsampx*dsampy);

dsampvecx = [ceil(dsampx/2):dsampx:lx];
dsampvecy = [ceil(dsampy/2):dsampy:ly];

E = zeros(lx,ly,nframes);

% define temporal filter
filtord = min(floor(log(.0001)/log(alpha)),100);
expfilt2 = alpha.^abs([-filtord:filtord]);
expfilt2 = expfilt2/sum(expfilt2);

saved_iterations = unique(find((mod([1:niter],iter_step)==0) | ([1:niter]==niter)));
sample_est_hist = zeros(lx,ly, length(saved_iterations));
sample_est = ones(lx,ly,nframes);

tic
ll=0;
for jiter = 1:niter
   %% deconSTORM iteration
   sample_est_down = cconv2_PSFk(sample_est, PSFk);
   sample_est_down = sample_est_down(dsampvecx,dsampvecy,:) + r;
   EDown = mov ./ sample_est_down;
   E(dsampvecx,dsampvecy,:) = EDown;
   Econv = cconv2_PSFk(E, PSFk);
   
   sample_est_curr = padarray(sample_est,[0,0,filtord],0,'pre');
   gam = filter(expfilt2, 1, sample_est_curr, [], 3);
   gam = gam(:,:,filtord+1:end);
   samplecurr = mean(sample_est,3);
   samplecurr = repmat(samplecurr,[1,1,nframes]);
   gam = min(gam + beta*samplecurr,1);
   gam = -log(gam) / gfactor;
   
   sample_est = sample_est .* Econv ./ (a+gam);
   
   %% Save intermediate iterations
   if (mod(jiter,iter_step)==0) || (jiter==niter)
      ll=ll+1;
      sample_est_hist(:,:,ll) = squeeze(mean(sample_est,3));
      tottime = toc;
      if verbose
         fprintf('deconSTORM: iteration %d/%d; time=%3.3g ms/iteration/frame; remaining time=%5.5g s.\n',...
            jiter,niter, 1000*tottime/jiter/nframes, tottime*(niter-jiter)/jiter)
      end
      
      %% Save output
      if ~isempty(fileout)
         sample_est_frames = sample_est;
         sample_est_mean = mean(sample_est,3);
         
         save(fileout,'sample_est_mean','sample_est_frames','sample_est_hist','saved_iterations');
         if verbose
            fprintf('  Saved current output to %s.\n',fileout)
         end
      end
   end
end

sample_est_hist = reshape(sample_est_hist, lx, ly, ceil(niter/iter_step));
sample_est_frames = reshape(sample_est, lx,ly,nframes);

% --------------- 
function op = cconv2_PSFk(x,PSFk)
%
% Convolution with periodic boundary conditions
%

if ndims(x)>=3
   xk = fft(fft(x,[],1),[],2);
   op = ifft(ifft(xk.*PSFk,[],1),[],2);
   op = fftshift(fftshift(real(op),1),2);
else
   op = fftshift(ifft2(fft2(x).*PSFk,'symmetric'));
end


