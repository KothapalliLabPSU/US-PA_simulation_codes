%% Overview
% This script simulates the back propagation of the PA affect generated
% acoustic wavefront

%% User inputs
Root = 'C:\Program Files\MATLAB\R2020a\toolbox\custom codes\custom codes';
n_wave = 1; % # of wavelengths
n_sources = 128;% # of sources, CONFIRM
wv_vect = [800];% wavelengths

%% Simulation
for i = 1:n_wave
    load(strcat('NIRdata',num2str(wv_vect(i)),'.mat'));
    
    figure;
    subplot(1,2,1);
    surf(NIRdata.x,NIRdata.y,NIRdata.ps2,'EdgeColor','none')
    NIRdata.y=NIRdata.y-30;
    NIR_size=size(NIRdata.x);
    
    PML_size = 40;          % size of the PML in grid points, can be changed accourding to user preferences
    Nx = NIR_size(2);  % number of grid points in the x (row) direction
    Ny = NIR_size(1);  % number of grid points in the y (column) direction
    dx = 0.1e-3;            % grid point spacing in the x direction [m]
    dy = 0.1e-3;            % grid point spacing in the y direction [m]
    kgrid = makeGrid(Nx, dx, Ny, dy);
    subplot(1,2,2);
    surf(kgrid.x,kgrid.y,zeros(size(kgrid.x)));
    A = interp2(NIRdata.x/1000,NIRdata.y/1000,NIRdata.ps2,kgrid.x,kgrid.y,'linear'); % Initial pressure input to k-wave - p0
    figure(3);
    subplot(1,2,2);
    h31=surf(kgrid.x,kgrid.y,A);
    set(h31, 'edgecolor','none')
    
    % define the properties of the propagation medium   
    medium.sound_speed = imrotate(sound_speed_mat,270);  % [m/s]
    medium.density = imrotate(density_mat,270);      % [kg/m^3]
    
    % Define initial pressure distribution in the domain
    source.p0=A;
    
    % smooth the initial pressure distribution and restore the magnitude
    source.p0 = smooth(kgrid, source.p0, true);
    
    % define a binary line sensor
    sensor.mask = zeros(Nx, Ny);
    sensor.mask = makeLine(Nx, Ny, [((Nx-1)/2-128) 1], [((Nx-1)/2+127) 1]);%changed from (-64:63) to (-128:127)
    sensor.mask(173:2:427,1)=0;
    
    % create the time array
    [kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);
    
    
    % set the input arguements: force the PML to be outside the computational
    % grid; switch off p0 smoothing within kspaceFirstOrder2D
    input_args = {'PMLInside', false, 'PMLSize', PML_size, 'PlotPML', false, 'Smooth', false};
    % run the simulation
     sensor_data_temp1 = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
     
    save(strcat('sensor_data',num2str(wv_vect(i)),'.mat'),'sensor_data_temp1');
end