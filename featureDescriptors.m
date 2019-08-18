%% CSC262: Feature Descriptors

%% Overview
% In this lab, we explore how feature descriptors can be found using a
% MOPS-like implementation. For simplicity, we leave out adjustments for
% multi-scale and orientation. We hope to learn how feature descriptors in
% an image can be effectively found without significant loss of
% information.


run ~/startup.m
I = im2double(imread('/home/weinman/courses/CSC262/images/kpdet1.tif'));

% Processing for Descriptors
detectI = kpdet(I); 
% Locate the rows and columns where keypoints were found 
[row, col] = find(detectI > 0); 
% Calculate the number of keypoints N
N = sum(detectI(:)); 
patches = zeros(N, 64); 
% Create the Gaussian kernel with scale of 5
sigma = 5; 
gauss = gkern(sigma^2);
% Apply the Gaussian kernel to blur input image
blurI = conv2(gauss, gauss, I, 'same'); 
% Downsample image by a factor of 5
downsampleFactor = 5;
downsampleI = blurI(1:downsampleFactor:end, 1:downsampleFactor:end); 


%% Extracting Patches
% In this section, we extract patches around our detected key points, which
% will later be normalized and used as our feature descriptors. 

% b1
k = N-1; % a value between 1 and N

% Find the scaled location of the keypoint by deviding by the downsample
% factor
downsampledRow = round(row(k)/downsampleFactor);
downsampledCol = round(col(k)/downsampleFactor);

% Find the coordinate of the upper-left point of the 8x8 patch 
sizeHalf = 4; 
upperleftRow = downsampledRow - sizeHalf; %approximate the upperleft coords
upperleftCol = downsampledCol - sizeHalf;

% Find the coordinate of the lower-right point of the 8x8 patch
lowerrightRow = downsampledRow + sizeHalf - 1;
lowerrightCol = downsampledCol + sizeHalf - 1;

% Extract patch from downsampled image
patch = downsampleI(upperleftCol:lowerrightCol,upperleftRow:lowerrightRow);

% b6
figure;
subplot(1,2,1); 
imagesc(patch,[0 1]); % Display the image with black 0 and white 1
colormap(gray);     % Render in grayscale
axis equal off;     % Use square pixels and turn off borders
title('Downsampled Patch');

% b7
% get upsampled version of patch
origupperleftRow = row(k) - sizeHalf*downsampleFactor;
origupperleftCol = col(k) - sizeHalf*downsampleFactor;
origlowerrightRow = row(k) + sizeHalf*downsampleFactor - 1;
origlowerrightCol = col(k) + sizeHalf*downsampleFactor - 1;

upsampledPatch = I(origupperleftCol:origlowerrightCol, ...
    origupperleftRow:origlowerrightRow);
subplot(1,2,2); 
imagesc(upsampledPatch,[0 1]); % Display the image with black 0 and white 1
colormap(gray);     % Render in grayscale
axis equal off;     % Use square pixels and turn off borders
title('Original Patch');

% B8
%%
% Above is the visualization of the down-sampled patch, along with the
% original patch for comparison. The down-sampled and original patches show
% changes in intensity in the similar locations, but the majority of the
% detail in the origianl patch is lost. In the down-sampled patch, the
% intensity changes are not as rapid and clear as in the original, due to
% blurring and down-sampling. However, general trends in brightness stay
% the same across both images; where there is a specific tone in the
% original patch, there are matching tones at similar locations in the
% down-sampled patch. Below is a reference for the key point location
% in our original image.


% b9
figure;
imshow(I);
hold on;
plot(row(k), col(k), '+r');
title('Reference of Key Point Location');

%% Patch Processing
% To avoid affine photometric variations, we need to normalize patches'
% intensities so that their mean is zero and their variance is one. 

% Normalize the bias by subtracting the mean of all pixels in the patch
mean8 = mean(patch(:)); 
norBiasPatch8 = patch - mean8; 

% Normalize the gain by dividing the bias-normalized version by the
% standard deviation
std8 = std(patch(:)); 
norGainPatch8 = norBiasPatch8/std8;

% Show 2 patches, normalized and not normalized
figure; 
subplot(1,2,1);
imagesc(patch, [0 1]); % Display the image with black 0 and white 1
colormap(gray);     % Render in grayscale
axis equal off;
title('Non-normalized patch'); 

subplot(1,2,2); 
imagesc(norGainPatch8); 
colormap(gray);     % Render in grayscale
axis equal off;
title('Normalized Patch');

%%
% Above is the visualization of our bias and gain normalized patch, along
% with the non-normalized patch. The normalized patch has greater contrast
% than the non-normalized patch. All the pixels in each patch have the same
% intensity relative to the other pixels in the patch, and the 2 patches
% are similar in structure. However, the range of brightnesses in the
% normalized patch expands outside of the 0-1 range, making the patch more 
% distinctive and easier to match.

% Processing Keypoints
descriptors = kpfeat(I, detectI); 
% Function to convert a 64x1 vector into a 8x8 matrix, in case a quick
% visualization of the patches in descriptors is required
% function [patch] = oneTo2(vector, size)
%     patch = zeros(size, size); 
%     for i = 1:size
%         patch(:, i) = vector((i-1)*size + 1: i*size); 
%     end
% end


%% Conclusion
% Through this lab, we investigated the process of finding, and normalizing
% feature descriptors, and then generalized it into a function. After
% finding the key points of an image, a feature descriptor patch can be
% extracted for matching by blurring and down-sampling an area around a key
% point to compensate for slight inaccuracies in locating key points.
% Lastly, the patch must be normalized to compensate for affine photometric
% variations between images. The generalization of this process will allow
% us to apply what we learned about feature descriptors to future image
% matching problems.

%% Acknowledgement
% The kpdet1 image was taken by Jerod Weinman in Curitiba, Brazil at the
% Pontifical Catholic University of Parana and are Copyright 2007, licensed
% under a Creative Commons Attribution-Noncommercial-Share Alike 4.0
% International License. Szeliski's textbook, "Computer Vision: Algorithms
% and Applications" was referred to for understanding the purpose of
% normalizing patches. The function gkern was provided by and used at the
% courtesy of Jerod Weinman. The kpdet function used in this lab was built
% in the previous lab by a group member. We also utilized directions,
% information, and code snippets from the Feature Descriptors Lab text
% written by Jerod Weinman for CSC 262: Computer Vision this semester.
% Finally, we also consulted lab write up guidelines published on the
% course website by Jerod Weinman.
