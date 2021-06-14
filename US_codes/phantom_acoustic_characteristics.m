%% Overview 
% This script is used to generate the acoustically hetergeneous speed and 
% density maps for a given in-silico phantom image.
%
% USER INPUT: phantom image as png file, line 10 
% OUTPUT: density and sound speed maps, to be used in simulation code next

clc; clear all; close all;
%% Loading the phantom 
phantom = rgb2gray(imread('example_phantom\phantom_image.png')); % include phantom image here
[rows, cols] = size(phantom);

figure; imshow(phantom,[]);

%% initializing sound speed and density maps
sound_speed_mat = zeros([rows cols]); % m/s
density_mat = zeros([rows cols]); % kg/m^3

%% texture dictionary
texture_dict = rgb2gray(imread('contrast_dictionary.png'));
figure; imshow(texture_dict,[]);

%% brightness levels 
% class 0 ... class 11 represent the texture brightnesses in ascending
% order, class 0 is darkest tissue texture, class 11 is brightest texture

% n_levels is the number of levels between the left and right boundaries of
% each contrast class, it can be varied, does not significantly change
% texture of each class, keep constant for each class
n_levels = 20;

class_0 = linspace(.99975, 1.00025, n_levels); % 0.05% variation, darkest texture
class_1 = linspace(.9995, 1.0005, n_levels); % 0.1% variation
class_2 = linspace(.999, 1.001, n_levels); % 0.2% variation
class_3 = linspace(.9985, 1.0015, n_levels); % 0.3% variation
class_4 = linspace(.998, 1.002, n_levels); % 0.4% variation
class_5 = linspace(.9975, 1.0025, n_levels); % 0.5% variation
class_6 = linspace(.997, 1.003, n_levels); % 0.6% variation
class_7 = linspace(.9965, 1.0035, n_levels); % 0.7% variation
class_8 = linspace(.996, 1.004, n_levels); % 0.8% variation
class_9 = linspace(.9955, 1.0045, n_levels); % 0.9% variation
class_10 = linspace(.995, 1.005, n_levels); % 1% variation
class_11 = linspace(.99, 1.01, n_levels); % 2% variation, brightest texture


%% assigning sound speed and density values
%
% define the sound speed and density of each tissue class
%
% assign a contrast level to each tissue class in the phantom (ex: pxl
% intensity 100 = class_1 contrast)
% randsample(class_1,1) <--- assign contrast classes in randsample function
%
% contrast levels can be either referenced from the contrast dictionary in the instruction
% manual, and also will be loaded at the beginning of this script
%
% below code randomly assigns values from the assigned texture vectors
% multiplied by tissue density to each pixel in a given tissue class.
% doing so creates density heterogeneity for the tissue class.
%
% default values below are for the example case
%
% user is free to redefine sound speed, density, contrast classes, and
% should add additional switch-cases for tissue classes if necessary
%

for i = 1:rows
    for j = 1:cols    
        switch phantom(i,j)
              
            % add additional cases for additional tissue classes as
            % necessary
            
            % adjust density and sound speed values as necessary
            
            case 100 % left target
                sound_speed_mat(i,j) = 1540;
                density = 1045;
                density_mat(i,j) = randsample(class_4,1) * density; % 0.4% density variation 
                
            case 150 % right target
                sound_speed_mat(i,j) = 1537;
                density = 1020;
                density_mat(i,j) = randsample(class_10,1) * density; % 1% density variation
             
            case 255 % background
                sound_speed_mat(i,j) = 1540;
                density = 1058; 
                density_mat(i,j) = randsample(class_11,1) * density; % 2% density variation
        end
    end
end
%% displaying density map and sound speed map

figure;imshow(density_mat, []);
figure;imshow(sound_speed_mat, []);

save("example_phantom\density_mat.mat", "density_mat");
save("example_phantom\sound_speed_mat.mat", "sound_speed_mat");

