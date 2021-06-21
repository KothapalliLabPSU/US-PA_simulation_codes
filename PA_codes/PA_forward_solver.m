clc;clear all;close all;

%% Overview 
% PA_forward_solver contains the finite element method (FEM) optical
% scattering simulation. using the previously defined optical absorbers, we
% can calculate the fluence properties of each. We can then calculate the
% resultant pressure/acoustic characteristics of each, based on the PA
% effect

%% User inputs 
Root = 'C:\Program Files\MATLAB\R2020a\toolbox\custom codes\custom codes\figures'; % adjust root to location of mesh_back and mesh_anom
wv_vect = [800]; % wavelength vector 
n_wave = 1; % # of wavelengths
n_sources = 500; % # of sources

%% Simulation
mesh_anom = load_mesh('mesh_anom');
mesh_back = load_mesh('mesh_back');

phantom_size = load('phantom_size.mat');
[rows cols] = size(phantom_size.phantom_size);

[mua, mus, kappa, E] = calc_mua_mus(mesh_anom,wv_vect); %calculates mua and mus for the mesh
[mua_back, mus_back, kappa_back, E_back] = calc_mua_mus(mesh_back,wv_vect); %same for back

forward_out = femdata_spectral(mesh_anom,100,wv_vect);%calculate the fluence
save_data(forward_out,'forward_out');

forward_out_back = femdata_spectral(mesh_back,100,wv_vect);%calculate the fluence for back. extra step
save_data(forward_out_back,'forward_out_back');

for i = 1:n_wave
    
    phi = full(abs(sum(forward_out.phi(:,(i-1)*n_sources+1:i*n_sources),2)));   % fluence calculation. for i =1 (:, 1 to 500)
                                                                                % for i = 2 --> (:, 501 to 1000) etc.
    mua_wv = mua(:,i);
    ps = mua_wv.*phi;
    
    phi_back = full(abs(sum(forward_out_back.phi(:,(i-1)*n_sources+1:i*n_sources),2))); % similarly for the back
    mua_wv_back = mua_back(:,i);
    ps_back = mua_wv_back.*phi_back;
    
    back_phi = plotim2_new(mesh_back,phi_back);
    fluence = plotim2_new(mesh_anom,phi);
    
    fluence2 = fluence;
    for j = 1:rows %fluence correction
        fluence2(:,j) = fluence(:,j).*(1./back_phi(:,(rows-1)/2));
    end
    figure;imshow(fluence2,[]);
    
    mua2 = plotim2_new(mesh_anom,mua_wv);
    ps2 = mua2.*fluence2;
    
    h0=figure;
    plotim2_new(mesh_anom,mua_wv);
    xt = get(gca, 'XTick');set(gca, 'FontName','Calibri','FontSize', 20, 'FontWeight', 'bold');
    title(strcat('mua_',num2str(wv_vect(i)),'nm'),'FontName','Calibri','FontSize',20,'FontWeight','bold'); box on;colorbar;axis equal;
    pathname = Root;    figfile2 = fullfile(pathname,strcat('mua_',num2str(wv_vect(i)),'.fig')); saveas(h0,figfile2);
    
    h1=figure;
    plotim2_new(mesh_anom,phi);
    xt = get(gca, 'XTick');set(gca, 'FontName','Calibri','FontSize', 20, 'FontWeight', 'bold');
    title(strcat('\phi_',num2str(wv_vect(i)),'nm'),'FontName','Calibri','FontSize',20,'FontWeight','bold'); box on;colorbar;axis equal;
    pathname = Root;    figfile2 = fullfile(pathname,strcat('\phi',num2str(wv_vect(i)),'.fig')); saveas(h1,figfile2);
    
    h2=figure;
    plotim2_new(mesh_anom,ps);
    xt = get(gca, 'XTick');set(gca, 'FontName','Calibri','FontSize', 20, 'FontWeight', 'bold');
    title(strcat('P_s: ',num2str(wv_vect(i)),'nm'),'FontName','Calibri','FontSize',20,'FontWeight','bold'); box on;colorbar;axis equal;
    pathname = Root;    figfile2 = fullfile(pathname,strcat('P_s_',num2str(wv_vect(i)),'.fig')); saveas(h2,figfile2);

    h3=figure;
    plotim2_new(mesh_back,abs(phi-phi_back));
    xt = get(gca, 'XTick');set(gca, 'FontName','Calibri','FontSize', 20, 'FontWeight', 'bold');
    title(strcat('\phi - \Phi b _',num2str(wv_vect(i)),'nm'),'FontName','Calibri','FontSize',20,'FontWeight','bold'); box on;colorbar;axis equal;
    pathname = Root;    figfile2 = fullfile(pathname,strcat('phi_phib_',num2str(wv_vect(i)),'.fig')); saveas(h3,figfile2);
    
    h4=figure;
    plotim2_new(mesh_back,abs(ps-ps_back));
    xt = get(gca, 'XTick');set(gca, 'FontName','Calibri','FontSize', 20, 'FontWeight', 'bold');
    title(strcat('P_s - P_b: ',num2str(wv_vect(i)),'nm'),'FontName','Calibri','FontSize',20,'FontWeight','bold'); box on;colorbar;axis equal;
    pathname = Root;    figfile2 = fullfile(pathname,strcat('P_s_Pb_',num2str(wv_vect(i)),'.fig')); saveas(h4,figfile2);

    % Export nodal and scattered pressure data
    NIRdata1 =zeros(length(mesh_anom.nodes),7);
    NIRdata1(:,1:3)=mesh_anom.nodes(:,1:3);
    NIRdata1(:,4)=phi;
    NIRdata1(:,5)=abs(phi-phi_back);
    NIRdata1(:,6)=ps;
    NIRdata1(:,7)=abs(ps-ps_back);
    exp1= sortrows(NIRdata1);
    NIRdata.x = reshape(exp1(:,1),rows,cols);
    NIRdata.y = reshape(exp1(:,2),rows,cols);
    NIRdata.z = reshape(exp1(:,3),rows,cols);
    NIRdata.phi = reshape(exp1(:,4),rows,cols);
    NIRdata.phidiff = reshape(exp1(:,5),rows,cols);
    NIRdata.ps = reshape(exp1(:,6),rows,cols);
    NIRdata.psdiff = reshape(exp1(:,7),rows,cols);
    NIRdata.ps2 = flipud(ps2);
    save(strcat('NIRdata',num2str(wv_vect(i)),'.mat'),'NIRdata');  % function form
end