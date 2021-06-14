clc;clear all;close all;

%% Overview 
% PA_forward_solver contains the finite element method (FEM) optical
% scattering simulation. using the previously defined optical absorbers, we
% can calculate the fluence properties of each. We can then calculate the
% resultant pressure/acoustic characteristics of each, based on the PA
% effect

%% User inputs 
Root = 'C:\Program Files\MATLAB\R2020a\toolbox\custom codes\custom codes\figures';
wv_vect = [800]; % wavelength vector 
n_wave = 1; % # of wavelengths
n_sources = 500; % # of sources

%% Simulation
mesh_anom = load_mesh('mesh_anom');
mesh_back = load_mesh('mesh_back');

[mua, mus, kappa, E] = calc_mua_mus(mesh_anom,wv_vect);
[mua_back, mus_back, kappa_back, E_back] = calc_mua_mus(mesh_back,wv_vect);

forward_out = femdata_spectral(mesh_anom,100,wv_vect); % calculate the fluence
save_data(forward_out,'forward_out');

forward_out_back = femdata_spectral(mesh_back,100,wv_vect); % calculates fluence for background mesh
save_data(forward_out_back,'forward_out_back');


for i = 1:n_wave
    
    phi = full(abs(sum(forward_out.phi(:,(i-1)*n_sources+1:i*n_sources),2))); % fluence calculation. for i =1 (:, 1 to 500)
                                                                                % for i = 2 --> (:, 501 to 1000) etc.  
    mua_wv = mua(:,i);
    ps = mua_wv.*phi; % getting the resultant pressure data from the PA effect
    
    phi_back = full(abs(sum(forward_out_back.phi(:,(i-1)*n_sources+1:i*n_sources),2)));
    mua_wv_back = mua_back(:,i);
    ps_back = mua_wv_back.*phi_back;  % similarly for the background (pressure data)
    
    back_phi = plotim2_new(mesh_back,phi_back);
    fluence = plotim2_new(mesh_anom,phi);
    
    fluence2 = fluence;
    for j = 1:601
        fluence2(:,j) = fluence(:,j).*(1./back_phi(:,300)); % fluence correction
    end
    figure;imshow(fluence2,[]);
    
    mua2 = plotim2_new(mesh_anom,mua_wv);
    ps2 = mua2.*fluence2; % corrected pressure data
    

    % Export nodal and scattered pressure data
    NIRdata1 =zeros(length(mesh_anom.nodes),7);
    NIRdata1(:,1:3)=mesh_anom.nodes(:,1:3);
    NIRdata1(:,4)=phi;
    NIRdata1(:,5)=abs(phi-phi_back);
    NIRdata1(:,6)=ps;
    NIRdata1(:,7)=abs(ps-ps_back);
    exp1= sortrows(NIRdata1);
    NIRdata.x = reshape(exp1(:,1),601,601);
    NIRdata.y = reshape(exp1(:,2),601,601);
    NIRdata.z = reshape(exp1(:,3),601,601);
    NIRdata.phi = reshape(exp1(:,4),601,601);
    NIRdata.phidiff = reshape(exp1(:,5),601,601);
    NIRdata.ps = reshape(exp1(:,6),601,601);
    NIRdata.psdiff = reshape(exp1(:,7),601,601);
    NIRdata.ps2 = flipud(ps2);
    save(strcat('NIRdata',num2str(wv_vect(i)),'.mat'),'NIRdata');  % function form
end