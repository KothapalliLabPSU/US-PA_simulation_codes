clc; clear;close all;
%% Overview
% This script takes the final pressure data from the detectors and forms
% the PA image

%% Image forming
data_str = load('sensor_data800.mat');
data_mat = squeeze(data_str.sensor_data_temp1);

% data_mat(:,1:200) = 0;

I0 = ones(1,5520).*10000;

d = linspace(0,6,5520);

alpha = 0.3; f = 10; 
I = I0.*exp(-0.115*alpha*f*d);
I_flip = fliplr(I);
figure;plot(I_flip);

data_comp = data_mat';

for i = 1:128
    data_comp(:,i) = data_comp(:,i).*I_flip';        
end

[rekon,rekonuncut] = rekon_OA_freqdom(fliplr(data_comp(1:4000,:)),100,.2,1.50,0,1,1,5,1);
rekon_reshape = imresize(rekon,[1024 1024]);
figure;imshow(JW_LogCompress(rekon_reshape,60),[]);colormap hot;