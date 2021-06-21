clc; clear all; close all;

%% Overview
% This script defines the background and foreground optical characteristics
% of the phantom
%% Loading mesh and phantom
mesh_spec = load_mesh('\mesh_1\mesh_1.rpt'); % USER input

phantom = rgb2gray(imread('phantom_image.png')); % USER input
[rows, cols] = size(phantom);

%% Background Characteristics 

map_x = zeros(1,cols+1);
map_y = zeros(1,rows+1);
map_x(1,:) = -(cols/20):0.1:(cols/20); % change division to floor division in the phantom if odd
map_y(1,:) = 0:0.1:(rows/10); % change division to floor division in the phantom if odd
map_y = fliplr(map_y);


mesh_back = mesh_spec;

for i = 1:length(mesh_spec.nodes)
    x_cord = mesh_spec.nodes(i,1);
    y_cord = mesh_spec.nodes(i,2);
    
    col_ind = find(abs(map_x-x_cord)<0.01);
    row_ind = find(abs(map_y-y_cord)<0.01);
    
    map_val = 32; % modify as necessary    
    
    if (map_val == 32) % background tissue, modify as necessary
        mesh_back.region(i,1) = 1; % regions are entirely optional, refer NIRFAST documentation for details
        mesh_back.sa(i,1) = 1.0; % scattering amplitude
        mesh_back.sp(i,1) = 1.0; % scattering power
        
        % for the concentration, pick which column the tissue is in based on the order of the tissues in the .excoef and .param files - ex: (soft tissue being the 3rd tissue listed, its concentration values will be stored in the 3rd column)
        mesh_back.conc(i,1) = 0; % concentration of HbO 
        mesh_back.conc(i,2) = 0; % concentration of prostate
        mesh_back.conc(i,3) = 1; % concentration of soft tissue
        mesh_back.conc(i,4) = 0; % concentration of bladder
        mesh_back.ri(i,1) = 1.33; % refractive index
        mesh_back.c(i,1) = 3e11/1.33; % speed of light through tissue 
    end   
        
end

plotmesh(mesh_back,'spec');
save_mesh(mesh_back,'mesh_back');

%% Foreground Characteristics

phantom_map = zeros([rows cols]);

if (mod(rows, 2) == 0) && (mod(cols, 2) == 0) % adjusting the phantom if evenly dimensioned (ex: 600x600)
    phantom_map = zeros([rows+1, cols+1]);
    phantom_map(1:rows,1:cols) = phantom;
    phantom_map(rows+1,1:cols) = phantom(rows,:);
    phantom_map(1:rows,cols+1) = phantom(:,cols);
    phantom_map(rows+1,cols+1) = phantom(rows,cols);
end

mesh_anom = mesh_spec;

for i = 1:length(mesh_spec.nodes)
    x_cord = mesh_spec.nodes(i,1);
    y_cord = mesh_spec.nodes(i,2);
    
    col_ind = find(abs(map_x-x_cord)<0.01);
    row_ind = find(abs(map_y-y_cord)<0.01);
    
    map_val = prostate_map(row_ind,col_ind);
    
    if map_val == 255 %HbO
        mesh_anom.region(i,1) = 1; % regions are entirely optional, refer NIRFAST documentation for details
        mesh_anom.sa(i,1) = 0.00001; % scattering amplitude 
        mesh_anom.sp(i,1) = 0.66; % scattering power
        
        % for the concentration, pick which column the tissue is in based on the order of the tissues in the .excoef and .param files - ex: (soft tissue being the 3rd tissue listed, its concentration values will be stored in the 3rd column)
        mesh_anom.conc(i,1) = 1;
        mesh_anom.conc(i,2) = 0;
        mesh_anom.conc(i,3) = 0;
        mesh_anom.conc(i,4) = 0;
        
        mesh_anom.ri(i,1) = 1.3; % refractive index
        mesh_anom.c(i,1) = 3e11/1.3; % speed of light thru tisue medium
    end
    
    if (map_val == 100) % Prostate
        mesh_anom.region(i,1) = 2;
        mesh_anom.sa(i,1) = 0.19;
        mesh_anom.sp(i,1) = 1;
        mesh_anom.conc(i,1) = 0;
        mesh_anom.conc(i,2) = 1;
        mesh_anom.conc(i,3) = 0;
        mesh_anom.conc(i,4) = 0;
        
        mesh_anom.ri(i,1) = 1.3;
        mesh_anom.c(i,1) = 3e11/1.3;
    end
    
    if (map_val == 32) %Tissue       
        mesh_anom.region(i,1) = 3;
        mesh_anom.sa(i,1) = 1.0;
        mesh_anom.sp(i,1) = 1.0;
        mesh_anom.conc(i,1) = 0;
        mesh_anom.conc(i,2) = 0;
        mesh_anom.conc(i,3) = 1;
        mesh_anom.conc(i,4) = 0;
        mesh_anom.ri(i,1) = 1.33;
        mesh_anom.c(i,1) = 3e11/1.33;
    end   
    
    if map_val == 150 %Bladder
        mesh_anom.region(i,1) = 4;
        mesh_anom.sa(i,1) = 1.6;
        mesh_anom.sp(i,1) = 1.0;
        mesh_anom.conc(i,1) = 0;
        mesh_anom.conc(i,2) = 0;
        mesh_anom.conc(i,3) = 0;
        mesh_anom.conc(i,4) = 1;
        mesh_anom.ri(i,1) = 1.2;
        mesh_anom.c(i,1) = 3e11/1.2;
    end     
end

plotmesh(mesh_anom,'spec');
save_mesh(mesh_anom,'mesh_anom');

abs_map = zeros(rows,cols);
abs_map(phantom==0) = 0.4;
abs_map(phantom==32) = 0.001;
abs_map(phantom==100) = 0.003;
abs_map(phantom==150) = 0.002;
figure;imshow(abs_map,[]);axis on;


