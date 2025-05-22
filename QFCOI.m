function [QFCOI_value] = QFCOI(ref_image, test_image)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Implementation of objective QFCOI method for compressed omnidirectional images

% If you use our method, please cite the following paper:
% M. Simka, L. Polak, A. Zizien and K. Fliegel, "QFCOI - Enhancing Objective
% Quality Assessment for Compressed Omnidirectional Images with Fusion of
% Measure,"

% Copyright (c) 2025 Marek Simka
% [marek.sim2@gmail.com]

% Permission to use, copy, modify, and/or distribute this software for any
% purpose with or without fee is hereby granted, provided that the above
% copyright notice and this permission notice appear in all copies.
%
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

% This software also uses the sources:
%
% Copyright (c) 2007 Visual Communications Lab, Cornell University
%
% Copyright(c) 2003 Zhou Wang
%  Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%  quality assessment: From error visibility to structural similarity"
%  IEEE Transactios on Image Processing, vol. 13, no. 4, pp.600-612,
%  Apr. 2004.
%
%  Lin Zhang, Lei Zhang, Xuanqin Mou, and David Zhang,"FSIM: a feature similarity
%  index for image qualtiy assessment", IEEE Transactions on Image Processing, 
%  vol. 20, no. 8, pp. 2378-2386, 2011.
%
% Copyright (c) 1996-2009 Peter Kovesi
%  Peter Kovesi, "Image Features From Phase Congruency". Videre: A
%  Journal of Computer Vision Research. MIT Press. Volume 1, Number 3,
%  Summer 1999 http://mitpress.mit.edu/e-journals/Videre/001/v13.html
%
% Copyright (c) 2005 The University of Texas at Austin
%  H. R. Sheikh and A. C. Bovik, "Image Information and Visual Quality"., IEEE Transactions on Image Processing, (to appear).
%  Download manuscript draft from http://live.ece.utexas.edu in the Publications link.
%
% Particular details are specified individually for each code in the script.

%%% Input:
% (1) ref_image - reference omnidirectional image uncompressed
% (2) test_image - distorted (compressed) omnidirectional image
%%% Usage:
% result = QFCOI(ref_image,test_image);
% see example code
 

%% Preprocessing  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Inspired by preprocess_metrix_mux.m for MS-SSIM and VIFp
%%%  Author:    Matthew Gaubatz
%%%  Version:   1.1
%%%  Date:      04/15/2007
%%%  Copyright (c) 2007 Visual Communications Lab, Cornell University
%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ref_image_preprocess = ref_image;
if ischar(ref_image_preprocess)
    try
        ref_image_preprocess = imread(ref_image_preprocess);
    catch
        error(sprintf('Error [preprocess]: Input ''%s'' is not a valid image or file.', ref_image_preprocess));
    end
end

if ~isa(ref_image_preprocess, 'double')
    ref_image_preprocess = double(ref_image_preprocess);
end

[H1, W1, D1] = size(ref_image_preprocess);

if D1 ~= 1
    if D1 == 3
        ref_image_preprocess = 0.29900 * ref_image_preprocess(:,:,1) + 0.58700 * ref_image_preprocess(:,:,2) + 0.11400 * ref_image_preprocess(:,:,3);
    else
        error('Error [preprocess]: ref_image_preprocess is not a recognized color format.');
    end
end

test_image_preprocess = test_image;
if ischar(test_image_preprocess)
    try
        test_image_preprocess = imread(test_image_preprocess);
    catch
        error(sprintf('Error: Input ''%s'' is not a valid image or file.', test_image_preprocess));
    end
end

if ~isa(test_image_preprocess, 'double')
    test_image_preprocess = double(test_image_preprocess);
end

[H2, W2, D2] = size(test_image_preprocess);

if D2 ~= 1
    if D2 == 3
        test_image_preprocess = 0.29900 * test_image_preprocess(:,:,1) + 0.58700 * test_image_preprocess(:,:,2) + 0.11400 * test_image_preprocess(:,:,3);
    else
        error('Error: test_image_preprocess is not a recognized color format.');
    end
end

%% QFCOI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSSSIM_value = mssim_index(ref_image_preprocess, test_image_preprocess); 
[~, FSIMc_value] = FeatureSIM(ref_image, test_image);
VIFp_value = vifp_mscale(ref_image_preprocess, test_image_preprocess); 
PIQE_value = piqe(test_image);
    
QFCOI_weights = [0.8913, 0.1029, 0.0053, 0.0005]; % weight coefficients for fusion

QFCOI_max = QFCOI_weights(1) * 1 + QFCOI_weights(2) * 1 + QFCOI_weights(3) * 1 + QFCOI_weights(4) * 100; % QFCOI_min = 0;
QFCOInonscaled = QFCOI_weights(1) * MSSSIM_value + QFCOI_weights(2) * FSIMc_value + QFCOI_weights(3) * VIFp_value + QFCOI_weights(4) * PIQE_value;
QFCOI_value = (QFCOInonscaled) / (QFCOI_max); % min-max normalization as QFCOI_value = (QFCOInonscaled - QFCOI_min) / (QFCOI_max - QFCOI_min);

end

%% MS-SSIM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mssim, comp, detail] = mssim_index(img1, img2, nlevs,K,lpf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  FUNCTION:  mssim_index
%%%
%%%  INPUTS:    img1                - original image data 
%%%
%%%             img2                - modified image data to be compared with
%%%                                   original image
%%%
%%%             nlevs               - number of levels for analysis 
%%%                                   (future parameter, currently only nlevs=5 
%%%                                   is supported)
%%%
%%%             K                   - see ssim_index
%%%
%%%             lpf                 - impulse response of low-pass filter to use
%%%                                   (defaults to 9/7 biorthogonal wavelet 
%%%                                   filters)
%%%
%%%  OUTPUTS:   mssim               - mssim value
%%%
%%%             comp                - vector containing mssim components as
%%%                                   [mean, contrast, cross-correlation]
%%%
%%%             detail              - vector containing mssim components at EACH
%%%                                   scale (same structure as comp)
%%%
%%%  CHANGES:   NONE
%%%             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Author:    David Rouse
%%%  Version:   1.0
%%%  Date:      01/05/2007
%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Copyright (c) 2007 Visual Communications Lab, Cornell University
%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Permission to use, copy, modify, distribute and sell this software
%%%  and its documentation for any purpose is hereby granted without fee,
%%%  provided that the above copyright notice appear in all copies and
%%%  that both that copyright notice and this permission notice appear
%%%  in supporting documentation.  VCL makes no representations about
%%%  the suitability of this software for any purpose.
%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  DISCLAIMER:
%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  The code provided hereunder is provided as is without warranty
%%%  of any kind, either express or implied, including but not limited
%%%  to the implied warranties of merchantability and fitness for a
%%%  particular purpose.  The author(s) shall in no event be liable for
%%%  any damages whatsoever including direct, indirect, incidental,
%%%  consequential, loss of business profits or special damages.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Multi-scale SSIM
% No need to specify first or second image as original and distorted
%  Input:
%   img1 - first image
%   img2 - second image
%
%  Output:
%   mssim - mssim value
%   comp - vector containing mssim components as
%           [mean, contrast, cross-correlation]
%   detail - vector containing mssim components at EACH scale (same
%   structure as comp)

nlevs = 5;

try
    K = K;
catch
    % Default SSIM parameters (assumes L=255)
    K = [0.01 0.03];
end

try
    lpf = lpf;
catch
    % Use Analysis Low Pass filter from Biorthogonal 9/7 Wavelet
    lod = [0.037828455507260; -0.023849465019560;  -0.110624404418440; ...
        0.377402855612830; 0.852698679008890;   0.377402855612830;  ...
        -0.110624404418440; -0.023849465019560; 0.037828455507260];
    lpf = lod*lod';
    lpf = lpf/sum(lpf(:));
end

img1 = double(img1);
img2 = double(img2);

window = fspecial('gaussian',11,1.5);
window = window/sum(sum(window));

ssim_v = zeros(nlevs,1);
ssim_r = zeros(nlevs,1);

% Scale 1 is the original image
[junk junk junk comp_ssim] = ssim_index_modified(img1,img2,K);
ssim_v(1) = comp_ssim(2);
ssim_r(1) = comp_ssim(3);

% Compute SSIM for scales 2 through 5
for s=1:nlevs-1    
    % Low Pass Filter
    img1 = imfilter(img1,lpf,'symmetric','same');
    img2 = imfilter(img2,lpf,'symmetric','same');
    
    img1 = img1(1:2:end,1:2:end);
    img2 = img2(1:2:end,1:2:end);

    [junk junk junk comp_ssim] = ssim_index_modified(img1,img2,K);
    ssim_m = comp_ssim(1);         % Mean Component only needed for scale 5
    ssim_v(s+1) = comp_ssim(2);
    ssim_r(s+1) = comp_ssim(3);   
end

alpha = 0.1333;
beta = [0.0448 0.2856 0.3001 0.2363 0.1333]';
gamma = [0.0448 0.2856 0.3001 0.2363 0.1333]';

comp = [ssim_m^alpha prod(ssim_v.^beta) prod(ssim_r.^gamma)];

detail = [ssim_m ssim_v' ssim_r'];

mssim = prod(comp);

return

end



function [ssim, ssim_map, map_obj, composite_mean_vec] = ssim_index_modified(img1, img2, K, window, L)
%%%MM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%MM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%MM%%% SSIM Index, Modified Version
%%%MM%%%
%%%MM%%% The authors of this modified version are with the Visual 
%%%MM%%% Communications Laboratory, in department of Electrical and 
%%%MM%%% Computer Engineering at Cornell University.  The original 
%%%MM%%% code was written by Zhou Wang at NYU. All of the information
%%%MM%%% and documentation provided with the original code is included
%%%MM%%% below.  For more information, see the following website:
%%%MM%%% 
%%%MM%%% http://www.ece.uwaterloo.ca/~z70wang/research/ssim/
%%%MM%%%
%%%MM%%% This version of the SSIM index computation routine is identical
%%%MM%%% to the original routine, except that it returns two additional
%%%MM%%% arguments:
%%%MM%%%
%%%MM%%% (1) 'map_obj', which has fields M, V and R which respectively
%%%MM%%%     represent the luminance, contrast and structure components
%%%MM%%%     used during the SSIM index computation, and
%%%MM%%%
%%%MM%%% (2) 'component_mean_vec', which contains the mean values of the
%%%MM%%%     of 'map_obj' fields.
%%%MM%%%
%%%MM%%% This routine is used to help compute the version of the MSSIM
%%%MM%%% index, a multi-scale version of the SSIM index, included with
%%%MM%%% this package.
%%%MM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%MM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%MM%%% Authors:   David Rouse and Matthew Gaubatz
%%%MM%%% Version:   1.0
%%%MM%%% Date:      01/05/2007
%%%MM%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%MM%%% Copyright (c) 2007 Visual Communications Lab, Cornell University 
%%%MM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%MM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%SSIM Index, Version 1.0
%Copyright(c) 2003 Zhou Wang
%All Rights Reserved.
%
%The author is with Howard Hughes Medical Institute, and Laboratory
%for Computational Vision at Center for Neural Science and Courant
%Institute of Mathematical Sciences, New York University.
%
%----------------------------------------------------------------------
%Permission to use, copy, or modify this software and its documentation
%for educational and research purposes only and without fee is hereby
%granted, provided that this copyright notice and the original authors'
%names appear on all copies and supporting documentation. This program
%shall not be used, rewritten, or adapted as the basis of a commercial
%software or hardware product without first obtaining permission of the
%authors. The authors make no representations about the suitability of
%this software for any purpose. It is provided "as is" without express
%or implied warranty.
%----------------------------------------------------------------------
%
%This is an implementation of the algorithm for calculating the
%Structural SIMilarity (SSIM) index between two images. Please refer
%to the following paper:
%
%Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%quality assessment: From error visibility to structural similarity"
%IEEE Transactios on Image Processing, vol. 13, no. 4, pp.600-612,
%Apr. 2004.
%
%Kindly report any suggestions or corrections to zhouwang@ieee.org
%
%----------------------------------------------------------------------
%
%Input : (1) img1: the first image being compared
%        (2) img2: the second image being compared
%        (3) K: constants in the SSIM index formula (see the above
%            reference). defualt value: K = [0.01 0.03]
%        (4) window: local window for statistics (see the above
%            reference). default widnow is Gaussian given by
%            window = fspecial('gaussian', 11, 1.5);
%        (5) L: dynamic range of the images. default: L = 255
%
%Output: (1) mssim: the mean SSIM index value between 2 images.
%            If one of the images being compared is regarded as 
%            perfect quality, then mssim can be considered as the
%            quality measure of the other image.
%            If img1 = img2, then mssim = 1.
%        (2) ssim_map: the SSIM index map of the test image. The map
%            has a smaller size than the input images. The actual size:
%            size(img1) - size(window) + 1.
%
%Default Usage:
%   Given 2 test images img1 and img2, whose dynamic range is 0-255
%
%   [mssim ssim_map] = ssim_index(img1, img2);
%
%Advanced Usage:
%   User defined parameters. For example
%
%   K = [0.05 0.05];
%   window = ones(8);
%   L = 100;
%   [mssim ssim_map] = ssim_index(img1, img2, K, window, L);
%
%See the results:
%
%   mssim                        %Gives the mssim value
%   imshow(max(0, ssim_map).^4)  %Shows the SSIM index map
%
%========================================================================


if (nargin < 2 | nargin > 5)
   ssim_index = -Inf;
   ssim_map = -Inf;
   return;
end

if (size(img1) ~= size(img2))
   ssim_index = -Inf;
   ssim_map = -Inf;
   return;
end

[M N] = size(img1);

if (nargin == 2)
   if ((M < 11) | (N < 11))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);	%
   K(1) = 0.01;								      % default settings
   K(2) = 0.03;								      %
   L = 255;                                  %
end

if (nargin == 3)
   if ((M < 11) | (N < 11))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);
   L = 255;
   if (length(K) == 2)
      if (K(1) < 0 | K(2) < 0)
		   ssim_index = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_index = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

if (nargin == 4)
   [H W] = size(window);
   if ((H*W) < 4 | (H > M) | (W > N))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   L = 255;
   if (length(K) == 2)
      if (K(1) < 0 | K(2) < 0)
		   ssim_index = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_index = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

if (nargin == 5)
   [H W] = size(window);
   if ((H*W) < 4 | (H > M) | (W > N))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   if (length(K) == 2)
      if (K(1) < 0 | K(2) < 0)
		   ssim_index = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_index = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
window = window/sum(sum(window));
img1 = double(img1);
img2 = double(img2);

mu1   = filter2(window, img1, 'valid');
mu2   = filter2(window, img2, 'valid');
mu1_sq = mu1.*mu1;
mu2_sq = mu2.*mu2;
mu1_mu2 = mu1.*mu2;
sigma1_sq = filter2(window, img1.*img1, 'valid') - mu1_sq;
sigma2_sq = filter2(window, img2.*img2, 'valid') - mu2_sq;
sigma1 = real(sqrt(sigma1_sq));
sigma2 = real(sqrt(sigma2_sq));
sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;

if (C1 > 0 && C2 > 0)
   M = (2*mu1_mu2 + C1)./(mu1_sq + mu2_sq + C1);
   V = (2*sigma1.*sigma2 + C2)./(sigma1_sq + sigma2_sq + C2);
   R = (sigma12 + C2/2)./(sigma1.*sigma2+C2/2);
   ssim_map = M.*V.*R;
else
   ssim_ln = 2*mu1_mu2;
   ssim_ld = mu1_sq + mu2_sq;
   M = ones(size(mu1));
   index_l = (ssim_ld>0);
   M(index_l) = ssim_ln(index_l)./ssim_ld(index_l);
   
   ssim_cn = 2*sigma1.*sigma2;
   ssim_cd = sigma1_sq + sigma2_sq;
   V = ones(size(mu1));
   index_c = (ssim_cd>0);
   V(index_c) = ssim_cn(index_c)./ssim_cd(index_c);
   
   ssim_sn = sigma12;
   ssim_sd = sigma1.*sigma2;
   R = ones(size(mu1));
   index1 = sigma1>0;
   index2 = sigma2>0;
   index_s = (index1.*index2>0);
   R(index_s) = ssim_sn(index_s)./ssim_sd(index_s);
   index_s = (index1.*not(index2)>0);
   R(index_s) = 0;
   index_s = (not(index1).*index2>0);
   R(index_s) = 0;
   
   ssim_map = M.*V.*R;
end

ssim = mean2(ssim_map);

%%%MM%%%
%%%MM%%% [BEGIN CHANGE]
%%%MM%%%

%%%MM%%%
%%%MM%%% create per-component maps
%%%MM%%%
map_obj.M = M;
map_obj.V = V;
map_obj.R = R;

%%%MM%%%
%%%MM%%% compute per-component map mean values
%%%MM%%%
composite_mean_vec = [mean2(M) mean2(V) mean2(R)];

%%%MM%%%
%%%MM%%% [END CHANGE]
%%%MM%%%

return
end


%% FSIMc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FSIM, FSIMc] = FeatureSIM(imageRef, imageDis)
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the algorithm for calculating the
% Feature SIMilarity (FSIM) index between two images.
%
% Please refer to the following paper
%
% Lin Zhang, Lei Zhang, Xuanqin Mou, and David Zhang,"FSIM: a feature similarity
% index for image qualtiy assessment", IEEE Transactions on Image Processing, vol. 20, no. 8, pp. 2378-2386, 2011.
% 
%----------------------------------------------------------------------
%
%Input : (1) imageRef: the first image being compared
%        (2) imageDis: the second image being compared
%
%Output: (1) FSIM: is the similarty score calculated using FSIM algorithm. FSIM
%	     only considers the luminance component of images. For colorful images, 
%            they will be converted to the grayscale at first.
%        (2) FSIMc: is the similarity score calculated using FSIMc algorithm. FSIMc
%            considers both the grayscale and the color information.
%Note: For grayscale images, the returned FSIM and FSIMc are the same.
%        
%-----------------------------------------------------------------------
%
%Usage:
%Given 2 test images img1 and img2. For gray-scale images, their dynamic range should be 0-255.
%For colorful images, the dynamic range of each color channel should be 0-255.
%
%[FSIM, FSIMc] = FeatureSIM(img1, img2);
%-----------------------------------------------------------------------

[rows, cols] = size(imageRef(:,:,1));
I1 = ones(rows, cols);
I2 = ones(rows, cols);
Q1 = ones(rows, cols);
Q2 = ones(rows, cols);

if ndims(imageRef) == 3 %images are colorful
    Y1 = 0.299 * double(imageRef(:,:,1)) + 0.587 * double(imageRef(:,:,2)) + 0.114 * double(imageRef(:,:,3));
    Y2 = 0.299 * double(imageDis(:,:,1)) + 0.587 * double(imageDis(:,:,2)) + 0.114 * double(imageDis(:,:,3));
    I1 = 0.596 * double(imageRef(:,:,1)) - 0.274 * double(imageRef(:,:,2)) - 0.322 * double(imageRef(:,:,3));
    I2 = 0.596 * double(imageDis(:,:,1)) - 0.274 * double(imageDis(:,:,2)) - 0.322 * double(imageDis(:,:,3));
    Q1 = 0.211 * double(imageRef(:,:,1)) - 0.523 * double(imageRef(:,:,2)) + 0.312 * double(imageRef(:,:,3));
    Q2 = 0.211 * double(imageDis(:,:,1)) - 0.523 * double(imageDis(:,:,2)) + 0.312 * double(imageDis(:,:,3));
else %images are grayscale
    Y1 = imageRef;
    Y2 = imageDis;
end

Y1 = double(Y1);
Y2 = double(Y2);
%%%%%%%%%%%%%%%%%%%%%%%%%
% Downsample the image
%%%%%%%%%%%%%%%%%%%%%%%%%
minDimension = min(rows,cols);
F = max(1,round(minDimension / 256));
aveKernel = fspecial('average',F);

aveI1 = conv2(I1, aveKernel,'same');
aveI2 = conv2(I2, aveKernel,'same');
I1 = aveI1(1:F:rows,1:F:cols);
I2 = aveI2(1:F:rows,1:F:cols);

aveQ1 = conv2(Q1, aveKernel,'same');
aveQ2 = conv2(Q2, aveKernel,'same');
Q1 = aveQ1(1:F:rows,1:F:cols);
Q2 = aveQ2(1:F:rows,1:F:cols);

aveY1 = conv2(Y1, aveKernel,'same');
aveY2 = conv2(Y2, aveKernel,'same');
Y1 = aveY1(1:F:rows,1:F:cols);
Y2 = aveY2(1:F:rows,1:F:cols);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the phase congruency maps
%%%%%%%%%%%%%%%%%%%%%%%%%
PC1 = phasecong2(Y1);
PC2 = phasecong2(Y2);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the gradient map
%%%%%%%%%%%%%%%%%%%%%%%%%
dx = [3 0 -3; 10 0 -10;  3  0 -3]/16;
dy = [3 10 3; 0  0   0; -3 -10 -3]/16;
IxY1 = conv2(Y1, dx, 'same');     
IyY1 = conv2(Y1, dy, 'same');    
gradientMap1 = sqrt(IxY1.^2 + IyY1.^2);

IxY2 = conv2(Y2, dx, 'same');     
IyY2 = conv2(Y2, dy, 'same');    
gradientMap2 = sqrt(IxY2.^2 + IyY2.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the FSIM
%%%%%%%%%%%%%%%%%%%%%%%%%
T1 = 0.85;  %fixed
T2 = 160; %fixed
PCSimMatrix = (2 * PC1 .* PC2 + T1) ./ (PC1.^2 + PC2.^2 + T1);
gradientSimMatrix = (2*gradientMap1.*gradientMap2 + T2) ./(gradientMap1.^2 + gradientMap2.^2 + T2);
PCm = max(PC1, PC2);
SimMatrix = gradientSimMatrix .* PCSimMatrix .* PCm;
FSIM = sum(sum(SimMatrix)) / sum(sum(PCm));

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the FSIMc
%%%%%%%%%%%%%%%%%%%%%%%%%
T3 = 200;
T4 = 200;
ISimMatrix = (2 * I1 .* I2 + T3) ./ (I1.^2 + I2.^2 + T3);
QSimMatrix = (2 * Q1 .* Q2 + T4) ./ (Q1.^2 + Q2.^2 + T4);

lambda = 0.03;

SimMatrixC = gradientSimMatrix .* PCSimMatrix .* real((ISimMatrix .* QSimMatrix) .^ lambda) .* PCm;
FSIMc = sum(sum(SimMatrixC)) / sum(sum(PCm));

return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ResultPC]=phasecong2(im)
% ========================================================================
% Copyright (c) 1996-2009 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby  granted, free of charge, to any  person obtaining a copy
% of this software and associated  documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% The software is provided "as is", without warranty of any kind.
% References:
%
%     Peter Kovesi, "Image Features From Phase Congruency". Videre: A
%     Journal of Computer Vision Research. MIT Press. Volume 1, Number 3,
%     Summer 1999 http://mitpress.mit.edu/e-journals/Videre/001/v13.html

nscale          = 4;     % Number of wavelet scales.    
norient         = 4;     % Number of filter orientations.
minWaveLength   = 6;     % Wavelength of smallest scale filter.    
mult            = 2;   % Scaling factor between successive filters.    
sigmaOnf        = 0.55;  % Ratio of the standard deviation of the
                             % Gaussian describing the log Gabor filter's
                             % transfer function in the frequency domain
                             % to the filter center frequency.    
dThetaOnSigma   = 1.2;   % Ratio of angular interval between filter orientations    
                             % and the standard deviation of the angular Gaussian
                             % function used to construct filters in the
                             % freq. plane.
k               = 2.0;   % No of standard deviations of the noise
                             % energy beyond the mean at which we set the
                             % noise threshold point. 
                             % below which phase congruency values get
                             % penalized.
epsilon         = .0001;                % Used to prevent division by zero.

thetaSigma = pi/norient/dThetaOnSigma;  % Calculate the standard deviation of the
                                        % angular Gaussian function used to
                                        % construct filters in the freq. plane.

[rows,cols] = size(im);
imagefft = fft2(im);              % Fourier transform of image

zero = zeros(rows,cols);
EO = cell(nscale, norient);       % Array of convolution results.                                 

estMeanE2n = [];
ifftFilterArray = cell(1,nscale); % Array of inverse FFTs of filters

% Pre-compute some stuff to speed up filter construction

% Set up X and Y matrices with ranges normalised to +/- 0.5
% The following code adjusts things appropriately for odd and even values
% of rows and columns.
if mod(cols,2)
    xrange = [-(cols-1)/2:(cols-1)/2]/(cols-1);
else
    xrange = [-cols/2:(cols/2-1)]/cols;	
end

if mod(rows,2)
    yrange = [-(rows-1)/2:(rows-1)/2]/(rows-1);
else
    yrange = [-rows/2:(rows/2-1)]/rows;	
end

[x,y] = meshgrid(xrange, yrange);

radius = sqrt(x.^2 + y.^2);       % Matrix values contain *normalised* radius from centre.
theta = atan2(-y,x);              % Matrix values contain polar angle.
                                  % (note -ve y is used to give +ve
                                  % anti-clockwise angles)
				  
radius = ifftshift(radius);       % Quadrant shift radius and theta so that filters
theta  = ifftshift(theta);        % are constructed with 0 frequency at the corners.
radius(1,1) = 1;                  % Get rid of the 0 radius value at the 0
                                  % frequency point (now at top-left corner)
                                  % so that taking the log of the radius will 
                                  % not cause trouble.

sintheta = sin(theta);
costheta = cos(theta);
clear x; clear y; clear theta;    % save a little memory

% Filters are constructed in terms of two components.
% 1) The radial component, which controls the frequency band that the filter
%    responds to
% 2) The angular component, which controls the orientation that the filter
%    responds to.
% The two components are multiplied together to construct the overall filter.

% Construct the radial filter components...

% First construct a low-pass filter that is as large as possible, yet falls
% away to zero at the boundaries.  All log Gabor filters are multiplied by
% this to ensure no extra frequencies at the 'corners' of the FFT are
% incorporated as this seems to upset the normalisation process when
% calculating phase congrunecy.
lp = lowpassfilter([rows,cols],.45,15);   % Radius .45, 'sharpness' 15

logGabor = cell(1,nscale);

for s = 1:nscale
    wavelength = minWaveLength*mult^(s-1);
    fo = 1.0/wavelength;                  % Centre frequency of filter.
    logGabor{s} = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));  
    logGabor{s} = logGabor{s}.*lp;        % Apply low-pass filter
    logGabor{s}(1,1) = 0;                 % Set the value at the 0 frequency point of the filter
                                          % back to zero (undo the radius fudge).
end

% Then construct the angular filter components...

spread = cell(1,norient);

for o = 1:norient
  angl = (o-1)*pi/norient;           % Filter angle.

  % For each point in the filter matrix calculate the angular distance from
  % the specified filter orientation.  To overcome the angular wrap-around
  % problem sine difference and cosine difference values are first computed
  % and then the atan2 function is used to determine angular distance.

  ds = sintheta * cos(angl) - costheta * sin(angl);    % Difference in sine.
  dc = costheta * cos(angl) + sintheta * sin(angl);    % Difference in cosine.
  dtheta = abs(atan2(ds,dc));                          % Absolute angular distance.
  spread{o} = exp((-dtheta.^2) / (2 * thetaSigma^2));  % Calculate the
                                                       % angular filter component.
end

% The main loop...
EnergyAll(rows,cols) = 0;
AnAll(rows,cols) = 0;

for o = 1:norient                    % For each orientation.
  sumE_ThisOrient   = zero;          % Initialize accumulator matrices.
  sumO_ThisOrient   = zero;       
  sumAn_ThisOrient  = zero;      
  Energy            = zero;      
  for s = 1:nscale,                  % For each scale.
    filter = logGabor{s} .* spread{o};   % Multiply radial and angular
                                         % components to get the filter. 
    ifftFilt = real(ifft2(filter))*sqrt(rows*cols);  % Note rescaling to match power
    ifftFilterArray{s} = ifftFilt;                   % record ifft2 of filter
    % Convolve image with even and odd filters returning the result in EO
    EO{s,o} = ifft2(imagefft .* filter);      

    An = abs(EO{s,o});                         % Amplitude of even & odd filter response.
    sumAn_ThisOrient = sumAn_ThisOrient + An;  % Sum of amplitude responses.
    sumE_ThisOrient = sumE_ThisOrient + real(EO{s,o}); % Sum of even filter convolution results.
    sumO_ThisOrient = sumO_ThisOrient + imag(EO{s,o}); % Sum of odd filter convolution results.
    if s==1                                 % Record mean squared filter value at smallest
      EM_n = sum(sum(filter.^2));           % scale. This is used for noise estimation.
      maxAn = An;                           % Record the maximum An over all scales.
    else
      maxAn = max(maxAn, An);
    end
  end                                       % ... and process the next scale

  % Get weighted mean filter response vector, this gives the weighted mean
  % phase angle.

  XEnergy = sqrt(sumE_ThisOrient.^2 + sumO_ThisOrient.^2) + epsilon;   
  MeanE = sumE_ThisOrient ./ XEnergy; 
  MeanO = sumO_ThisOrient ./ XEnergy; 

  % Now calculate An(cos(phase_deviation) - | sin(phase_deviation)) | by
  % using dot and cross products between the weighted mean filter response
  % vector and the individual filter response vectors at each scale.  This
  % quantity is phase congruency multiplied by An, which we call energy.

  for s = 1:nscale,       
      E = real(EO{s,o}); O = imag(EO{s,o});    % Extract even and odd
                                               % convolution results.
      Energy = Energy + E.*MeanE + O.*MeanO - abs(E.*MeanO - O.*MeanE);
  end

  % Compensate for noise
  % We estimate the noise power from the energy squared response at the
  % smallest scale.  If the noise is Gaussian the energy squared will have a
  % Chi-squared 2DOF pdf.  We calculate the median energy squared response
  % as this is a robust statistic.  From this we estimate the mean.
  % The estimate of noise power is obtained by dividing the mean squared
  % energy value by the mean squared filter value

  medianE2n = median(reshape(abs(EO{1,o}).^2,1,rows*cols));
  meanE2n = -medianE2n/log(0.5);
  estMeanE2n(o) = meanE2n;

  noisePower = meanE2n/EM_n;                       % Estimate of noise power.

  % Now estimate the total energy^2 due to noise
  % Estimate for sum(An^2) + sum(Ai.*Aj.*(cphi.*cphj + sphi.*sphj))

  EstSumAn2 = zero;
  for s = 1:nscale
    EstSumAn2 = EstSumAn2 + ifftFilterArray{s}.^2;
  end

  EstSumAiAj = zero;
  for si = 1:(nscale-1)
    for sj = (si+1):nscale
      EstSumAiAj = EstSumAiAj + ifftFilterArray{si}.*ifftFilterArray{sj};
    end
  end
  sumEstSumAn2 = sum(sum(EstSumAn2));
  sumEstSumAiAj = sum(sum(EstSumAiAj));

  EstNoiseEnergy2 = 2*noisePower*sumEstSumAn2 + 4*noisePower*sumEstSumAiAj;

  tau = sqrt(EstNoiseEnergy2/2);                     % Rayleigh parameter
  EstNoiseEnergy = tau*sqrt(pi/2);                   % Expected value of noise energy
  EstNoiseEnergySigma = sqrt( (2-pi/2)*tau^2 );

  T =  EstNoiseEnergy + k*EstNoiseEnergySigma;       % Noise threshold

  % The estimated noise effect calculated above is only valid for the PC_1 measure. 
  % The PC_2 measure does not lend itself readily to the same analysis.  However
  % empirically it seems that the noise effect is overestimated roughly by a factor 
  % of 1.7 for the filter parameters used here.

  T = T/1.7;        % Empirical rescaling of the estimated noise effect to 
                    % suit the PC_2 phase congruency measure
  Energy = max(Energy - T, zero);          % Apply noise threshold

  EnergyAll = EnergyAll + Energy;
  AnAll = AnAll + sumAn_ThisOrient;
end  % For each orientation
ResultPC = EnergyAll ./ AnAll;
return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOWPASSFILTER - Constructs a low-pass butterworth filter.
%
% usage: f = lowpassfilter(sze, cutoff, n)
% 
% where: sze    is a two element vector specifying the size of filter 
%               to construct [rows cols].
%        cutoff is the cutoff frequency of the filter 0 - 0.5
%        n      is the order of the filter, the higher n is the sharper
%               the transition is. (n must be an integer >= 1).
%               Note that n is doubled so that it is always an even integer.
%
%                      1
%      f =    --------------------
%                              2n
%              1.0 + (w/cutoff)
%
% The frequency origin of the returned filter is at the corners.
%
% See also: HIGHPASSFILTER, HIGHBOOSTFILTER, BANDPASSFILTER
%

% Copyright (c) 1999 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% October 1999
% August  2005 - Fixed up frequency ranges for odd and even sized filters
%                (previous code was a bit approximate)

function f = lowpassfilter(sze, cutoff, n)
    
    if cutoff < 0 || cutoff > 0.5
	error('cutoff frequency must be between 0 and 0.5');
    end
    
    if rem(n,1) ~= 0 || n < 1
	error('n must be an integer >= 1');
    end

    if length(sze) == 1
	rows = sze; cols = sze;
    else
	rows = sze(1); cols = sze(2);
    end

    % Set up X and Y matrices with ranges normalised to +/- 0.5
    % The following code adjusts things appropriately for odd and even values
    % of rows and columns.
    if mod(cols,2)
	xrange = [-(cols-1)/2:(cols-1)/2]/(cols-1);
    else
	xrange = [-cols/2:(cols/2-1)]/cols;	
    end

    if mod(rows,2)
	yrange = [-(rows-1)/2:(rows-1)/2]/(rows-1);
    else
	yrange = [-rows/2:(rows/2-1)]/rows;	
    end
    
    [x,y] = meshgrid(xrange, yrange);
    radius = sqrt(x.^2 + y.^2);        % A matrix with every pixel = radius relative to centre.
    f = ifftshift( 1 ./ (1.0 + (radius ./ cutoff).^(2*n)) );   % The filter
    return;

end

%% VIFp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function vifp=vifp_mscale(ref,dist)
% -----------COPYRIGHT NOTICE STARTS WITH THIS LINE------------
% Copyright (c) 2005 The University of Texas at Austin
% All rights reserved.
% 
% Permission is hereby granted, without written agreement and without license or royalty fees, to use, copy, 
% modify, and distribute this code (the source files) and its documentation for
% any purpose, provided that the copyright notice in its entirety appear in all copies of this code, and the 
% original source of this code, Laboratory for Image and Video Engineering (LIVE, http://live.ece.utexas.edu)
% at the University of Texas at Austin (UT Austin, 
% http://www.utexas.edu), is acknowledged in any publication that reports research using this code. The research
% is to be cited in the bibliography as:
% 
% H. R. Sheikh and A. C. Bovik, "Image Information and Visual Quality", IEEE Transactions on 
% Image Processing, (to appear).
% 
% IN NO EVENT SHALL THE UNIVERSITY OF TEXAS AT AUSTIN BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, 
% OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF THIS DATABASE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF TEXAS
% AT AUSTIN HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% THE UNIVERSITY OF TEXAS AT AUSTIN SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE DATABASE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS,
% AND THE UNIVERSITY OF TEXAS AT AUSTIN HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
% 
% -----------COPYRIGHT NOTICE ENDS WITH THIS LINE------------
% 
% This software release consists of a MULTISCALE PIXEL DOMAIN, SCALAR GSM implementation of the algorithm described in the paper:
% 
% H. R. Sheikh and A. C. Bovik, "Image Information and Visual Quality"., IEEE Transactions on Image Processing, (to appear).
% Download manuscript draft from http://live.ece.utexas.edu in the Publications link.
% 
% THE PIXEL DOMAIN ALGORITHM IS NOT DESCRIBED IN THE PAPER. THIS IS A COMPUTATIONALLY SIMPLER
% DERIVATIVE OF THE ALGORITHM PRESENTED IN THE PAPER
% 
% Input : (1) img1: The reference image as a matrix
%         (2) img2: The distorted image (order is important)
% 
% Output: (1) VIF the visual information fidelity measure between the two images
% 
% Default Usage:
%    Given 2 test images img1 and img2, whose dynamic range is 0-255
% 
%    vif = vifvec(img1, img2);
% 
% Advanced Usage:
%    Users may want to modify the parameters in the code. 
%    (1) Modify sigma_nsq to find tune for your image dataset.
% Email comments and bug reports to hamid.sheikh@ieee.org


sigma_nsq=2;

num=0;
den=0;
for scale=1:4
   
    N=2^(4-scale+1)+1;
    win=fspecial('gaussian',N,N/5);
    
    if (scale >1)
        ref=filter2(win,ref,'valid');
        dist=filter2(win,dist,'valid');
        ref=ref(1:2:end,1:2:end);
        dist=dist(1:2:end,1:2:end);
    end
    
    mu1   = filter2(win, ref, 'valid');
    mu2   = filter2(win, dist, 'valid');
    mu1_sq = mu1.*mu1;
    mu2_sq = mu2.*mu2;
    mu1_mu2 = mu1.*mu2;
    sigma1_sq = filter2(win, ref.*ref, 'valid') - mu1_sq;
    sigma2_sq = filter2(win, dist.*dist, 'valid') - mu2_sq;
    sigma12 = filter2(win, ref.*dist, 'valid') - mu1_mu2;
    
    sigma1_sq(sigma1_sq<0)=0;
    sigma2_sq(sigma2_sq<0)=0;
    
    g=sigma12./(sigma1_sq+1e-10);
    sv_sq=sigma2_sq-g.*sigma12;
    
    g(sigma1_sq<1e-10)=0;
    sv_sq(sigma1_sq<1e-10)=sigma2_sq(sigma1_sq<1e-10);
    sigma1_sq(sigma1_sq<1e-10)=0;
    
    g(sigma2_sq<1e-10)=0;
    sv_sq(sigma2_sq<1e-10)=0;
    
    sv_sq(g<0)=sigma2_sq(g<0);
    g(g<0)=0;
    sv_sq(sv_sq<=1e-10)=1e-10;
    
    
     num=num+sum(sum(log10(1+g.^2.*sigma1_sq./(sv_sq+sigma_nsq))));
     den=den+sum(sum(log10(1+sigma1_sq./sigma_nsq)));
    
end
vifp=num/den;

 end



