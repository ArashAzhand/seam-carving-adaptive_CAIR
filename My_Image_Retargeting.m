clc; clear; close all;

image_name = 'Dolls'; 
percentage_reduction = 0.5;

I = im2double(imread(['CV_Project\Samples dataset\' image_name '\' image_name '.png']));
DMap = im2double(imread(['CV_Project\Samples dataset\' image_name '\' image_name '_DMap.png']));
SMap = im2double(imread(['CV_Project\Samples dataset\' image_name '\' image_name '_SMap.png']));
GMap = computeEnergyMap(I);


alpha = calculate_adaptive_alpha(DMap); 

EMap = medfilt2((alpha) .* normalize(DMap) + (1-alpha) .* normalize(SMap)  , [5 5]);

importance_T = 0.75;
if(max(DMap(:)) > importance_T)
    SC_percentage = 0.3;
    res = seamCarving(I, DMap, SMap, SC_percentage * percentage_reduction);
    output_width = (1 - percentage_reduction) *  size(I , 2);
    res = imresize(res , [size(I , 1) , output_width] , 'bicubic');
    imwrite(res, strcat('CV_Project/Output/', image_name, '/', image_name, '_', num2str(percentage_reduction), '.png'));
else
    SC_percentage = 1;
    res = seamCarving(I, DMap, SMap, SC_percentage * percentage_reduction);
    imwrite(res, strcat('CV_Project/Output/', image_name, '/', image_name, '_', num2str(percentage_reduction), '.png'));
end

% Display the result
figure ,imshow([res, I], []);

function output_image = seamCarving(input_image, DMap , SMap,  SC_percentage)
   
    input_image = im2double(input_image);
    
    
    [height, width, channels] = size(input_image);
    
    output_width = width - (width * SC_percentage);
    num_vertical_seams = round(width - output_width);

    
    % Reduce width
    for i = 1:num_vertical_seams

        gray_image = im2double(rgb2gray(input_image));
    
        med_input = medfilt2(gray_image);
        G_emap = double(imgradient(med_input));
        S_emap = double(edge(med_input , 'sobel'));
        C_emap = double(edge(med_input , 'canny'));
    
        energy_map = 2 * normalize(S_emap) +  2 * normalize(C_emap) + 2 * normalize(G_emap) + ...
            3 * normalize(DMap) + 1 * normalize(SMap);
        % GMap = im2double(computeEnergyMap(input_image));
        % beta = calculate_adaptive_alpha(GMap);
        % energy_map = medfilt2(normalize(GMap) * (beta) + (1 - beta) * normalize(EMap) , [5,5]);
        
        seam = findVerticalSeam(energy_map);
        for j = 1:height
            input_image(j, seam(j), :) = [1, 0, 0]; % Red color
        end
        
        imshowpair(input_image,energy_map,"montage");
        % pause(0.05);

        [DMap , SMap ,energy_map, input_image] = removeVerticalSeam(input_image, seam, energy_map , DMap, SMap );
        
    end

    output_image = input_image;
end

function y = normalize(x)
%NORMALIZE Normalize the input array using min-max normalization
y = (x - min(x(:))) / (max(x(:)) - min(x(:)));
end

function energy_map = computeEnergyMap(input_image , EMap)
    gray_image = im2double(rgb2gray(input_image));
    
    med_input = medfilt2(gray_image);
    G_emap= imgradient(med_input);
    S_emap = edge(med_input , 'sobel');
    C_emap = edge(med_input , 'canny');

    energy_map = 0.3 * S_emap + 0.3 * C_emap + 0.4 * G_emap;

end

function seam = findVerticalSeam(energy_map)
    [height, width] = size(energy_map);
    M = energy_map;
    backtrack = zeros(height, width);

    for i = 2:height
        for j = 1:width
            if j == 1
                [val, idx] = min([M(i-1, j), M(i-1, j+1)]);
                backtrack(i, j) = idx + j - 1;
            elseif j == width
                [val, idx] = min([M(i-1, j-1), M(i-1, j)]);
                backtrack(i, j) = idx + j - 2;
            else
                [val, idx] = min([M(i-1, j-1), M(i-1, j), M(i-1, j+1)]);
                backtrack(i, j) = idx + j - 2;
            end
            M(i, j) = M(i, j) + val;
        end
    end

    % Backtrack to find the seam
    seam = zeros(height, 1);
    [~, seam(height)] = min(M(height, :));
    for i = height-1:-1:1
        seam(i) = backtrack(i+1, seam(i+1));
    end
end

function [dmap , smap, energy_map, input_image] = removeVerticalSeam(input_image, seam, energy_map , dmap , smap )
    [height, width, channels] = size(input_image);
    for i = 1:height
        for c = 1:channels
            input_image(i, seam(i):width-1, c) = input_image(i, seam(i)+1:width, c);
        end
        energy_map(i, seam(i):width-1) = energy_map(i, seam(i)+1:width);
        dmap(i, seam(i):width-1) = dmap(i, seam(i)+1:width);
        smap(i, seam(i):width-1) = smap(i, seam(i)+1:width);
    end
    input_image = input_image(:, 1:width-1, :);
    energy_map = energy_map(:, 1:width-1, :);
    dmap = dmap(:, 1:width-1, :);
    smap = smap(:, 1:width-1, :);
end

function alpha = calculate_adaptive_alpha(depth_map)
    % Convert depth_map to double if it's not already
    depth_map = im2double(depth_map);
    
    % Number of quantization levels
    K = round(max(depth_map(:)) * 256);
    
    % Quantize the depth map
    quantized_depth = round(depth_map * (K-1));
    
    % Calculate histogram
    [counts, bins] = histcounts(quantized_depth, K);
    
    % Remove bins with very few pixels (e.g., less than 0.1% of total pixels)
    threshold = numel(depth_map) * 0.001;
    valid_bins = counts > threshold;
    valid_counts = counts(valid_bins);
    valid_bins = bins(valid_bins);
    
    % Calculate mean
    mean_value = sum(valid_bins .* valid_counts) / sum(valid_counts);
    
    % Calculate dispersion (variance)
    D = sum(valid_counts .* (valid_bins - mean_value).^2);
    
    % Calculate D_max (assuming it's for a half black, half white image)
    D_max = sum(valid_counts .* ((K-1)/2).^2);
    
    % Calculate ratio
    ratio = D / D_max;
    
    % Set alpha range
    alpha_min = 0.2;
    alpha_max = 0.8;
    
    % Calculate alpha
    alpha = alpha_min + (alpha_max - alpha_min) * ratio;
    
    % Ensure alpha is within [alpha_min, alpha_max]
    alpha = max(alpha_min, min(alpha_max, alpha));
end