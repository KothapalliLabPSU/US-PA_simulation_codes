%% Overview
% This script simulates the back propagation of the PA affect generated
% acoustic wavefront

%% User inputs
Root = 'C:\Program Files\MATLAB\R2020a\toolbox\custom codes\custom codes'; % adjust root to location of sound and density matrices
n_wave = 1; % # of wavelengths
n_sources = 128;% # of sources, CONFIRM
wv_vect = [800];% wavelengths

%% Simulation
density_mat = load('density_mat.mat');
density_mat = density_mat.density_mat;

sound_speed_mat = load('sound_speed_mat.mat');
sound_speed_mat = sound_speed_mat.sound_speed_mat;

for i = 1:n_wave
    load(strcat('NIRdata',num2str(wv_vect(i)),'.mat'));%loading initial pressure generated from nirfast simulations (forward solver)
    
    figure;
    subplot(1,2,1);
    surf(NIRdata.x,NIRdata.y,NIRdata.ps2,'EdgeColor','none')%plotting ps at x , y to visualize
    NIRdata.y=NIRdata.y-30; %since we want it to move from -30 to 30
    NIR_size=size(NIRdata.x);

    PML_size = 150;          
    Nx = NIR_size(2);  % number of grid points in the x (row) direction (no. of pixels in the x direction)
    Ny = NIR_size(1);  % number of grid points in the y (column) direction (" " for y)
    dx = 0.1e-3;            % grid point spacing in the x direction [m]
    dy = 0.1e-3;            % grid point spacing in the y direction [m]
    kgrid = makeGrid(Nx, dx, Ny, dy); %defining the grid in k-wave, grid is same as mesh. 
    
    A = interp2(NIRdata.x/1000,NIRdata.y/1000,NIRdata.ps2,kgrid.x,kgrid.y,'linear'); % Initial pressure input to k-wave - p0, initial pressure mapped to k-grid
    figure(3);

    h31=surf(kgrid.x,kgrid.y,A);
    set(h31, 'edgecolor','none')
    % define the properties of the propagation medium   
    medium.sound_speed = imrotate(sound_speed_mat,270);  % [m/s] imrotate used to move bottom side of mesh to the left in grid
    medium.density = imrotate(density_mat,270);      % [kg/m^3]
    % Define initial pressure distribution in the domain
    source.p0=A; %stores the initial pressures in A to source.p0 source is the structure
    % smooth the initial pressure distribution and restore the magnitude
    source.p0 = smooth(kgrid, source.p0, true);
   
    % define a binary line sensor
    sensor.mask = zeros(Nx, Ny); 
    sensor.mask = makeLine(Nx, Ny, [((Nx-1)/2-128) 1], [((Nx-1)/2+127) 1]);%changed from (-64:63) to (-128:127)defines where the transducers are.
    sensor.mask(173:2:427,1)=0;

    [kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed); 
    input_args = {'PMLInside', false, 'PMLSize', PML_size, 'PlotPML', false, 'Smooth', false, 'RecordMovie', true, 'MovieName', 'wave_propagation', 'MovieArgs', {'FrameRate', 10}};

    sensor_data_temp1 = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:},'DataCast','gpuArray-single');% mimics propagation and records the pressures at the transducers defined in ln80
    sensor_data_temp1 = gather(sensor_data_temp1);
    save(strcat('sensor_data',num2str(wv_vect(i)),'.mat'),'sensor_data_temp1');
    sensor.time_reversal_boundary_data = sensor_data_temp1;
    PA_Image = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:},'DataCast','gpuArray-single');
    PA_Image = gather(PA_Image);
    log_image = JW_LogCompress(PA_Image, 40);
    figure; imshow(imrotate(PA_Image, 90), []);
    figure; imshow(imrotate(log_image, 90), []);
end