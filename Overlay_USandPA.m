% OVERVIEW: Code for overlay of PA image over US image
%
% Inputs: US and PA images in .png format
%
% Output: saves the overlaid USPA image in current working directory

clc;close all;clear;

US = imread('US.png'); % Reads the US image
PA = imread('PA.png'); % Reads the PA image

% Displaying US and PA images in gray and hot scales 
figure;imagesc(US);colormap gray;axis off;
figure;imagesc(PA);colormap hot;axis off;

% Setting up the overlay code
Threshold_PA = 0.05; % Threshold value for PA map to be displayed over US, in order to cut some
% noise in PA image
Overlay_perc = 0.75; % This parameter defines the percentage of transparency; 
% 1 = no US visible; 0 = no PA visible; intermediate values control
% transparency
alpha_mat = mat2gray(rgb2gray(PA));
alpha_mat(alpha_mat < Threshold_PA) = 0;
alpha_mat(alpha_mat > Threshold_PA) = Overlay_perc;

figure;
US2 = image(US);
axis off;
hold on;
PA2 = image(PA);
set(PA2,'AlphaData',alpha_mat);
saveas(gcf,'USPA_image.png'); %save the overlaid image in to current directory