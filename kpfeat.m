function [descriptors] = kpfeat(I,detectI)
% KPFEAT - finds feature descriptors of an interest point from an image.
%   I - The image to find the feature descriptors for.
%   detectI - a matrix containing the direction of the gradient of the
%             image at each point.
% returns a matrix containing the descriptors for each feature found in I.

    % Locate the rows and columns where keypoints were found 
    [row, col] = find(detectI > 0); 
    
    % Calculate the number of keypoints N
    N = sum(detectI(:)); 
    
    % Create the Gaussian kernel with scale of 5
    sigma = 5; 
    gauss = gkern(sigma^2);
    % Apply the Gaussian kernel to blur input image
    blurI = conv2(gauss, gauss, I, 'same'); 
    % Downsample image by a factor
    downsampleFactor = 5;
    downsampleI = blurI(1:downsampleFactor:end, 1:downsampleFactor:end);
    % Allocate space for the descriptor
    descriptors = zeros(N, 64);
    [hor, ver] = size(I); 
    % loop over all the keypoints in image to extract and normalize patches
    for i = 1:N
        % Find the scaled location of the keypoint by deviding by the
        % downsample factor
        downsampledRow = round(row(i)/downsampleFactor);
        downsampledCol = round(col(i)/downsampleFactor);
        % Find the coordinate of the upper-left point of the 8x8 patch 
        sizeHalf = 4; 
        upperleftRow = downsampledRow - sizeHalf; 
        upperleftCol = downsampledCol - sizeHalf;
        % Find the coordinate of the lower-right point of the 8x8 patch
        lowerrightRow = downsampledRow + sizeHalf - 1; 
        lowerrightCol = downsampledCol + sizeHalf - 1;
        % Check if the patch is out of bound of the image
        if (upperleftRow < 1 || upperleftCol < 1 ...
            || lowerrightRow > ver || lowerrightCol > hor)
            descriptors(i, :) = NaN;
            continue; 
        end
        
        % Extract patch from downsampled image
        patch = downsampleI(upperleftCol:lowerrightCol,upperleftRow:lowerrightRow);
        % Normalize the bias by subtracting the mean of all pixels in the
        % patch
        mean8 = mean(patch(:)); 
        norBiasPatch8 = patch - mean8; 
        % Normalize the gain by dividing the bias-normalized version by the
        % standard deviation
        std8 = std(patch(:)); 
        norGainPatch8 = norBiasPatch8/std8;
        % Setting descriptor for a feature point
        descriptors(i, :) = norGainPatch8(:);
    end
end

