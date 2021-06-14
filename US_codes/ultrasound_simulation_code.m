% OVERVIEW
% This script is used for simulating ultrasound using the k-Wave acoustics
% toolbox. 
%
% Inputs to this script are a density map and sound speed map of an
% in-silico phantom. Both maps will be automatically loaded in to the
% script.
%
% The primary output is the simulated ultrasound image of the in-silico
% phantom.

clc; clear all; close all;

density_mat = load('example_phantom\density_mat.mat'); % change path as necessary
density_map = flipud(density_mat.density_mat);

sound_speed_mat = load('example_phantom\sound_speed_mat.mat'); % change path as necessary
sound_speed_map = flipud(sound_speed_mat.sound_speed_mat);

[rows, cols] = size(density_map);

% simulation settings
DATA_CAST       = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
RUN_SIMULATION  = true;         % set to false to reload previous results instead of running simulation

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
pml_x_size = 5;                % [grid points]
pml_y_size = 5;                % [grid points]
pml_z_size = 0;                % [grid points]

% set total number of grid points not including the PML
sc = 1;
Nx = (cols+2*pml_x_size)/sc - 2*pml_x_size;     % [grid points]
Ny = (rows+2*pml_y_size)/sc - 2*pml_y_size;     % [grid points]
Nz = 4/sc - 2*pml_z_size;     % [grid points]

% set desired grid size in the x-direction not including the PML
x = 60e-3;                      % [m]

% calculate the spacing between the grid points
dx = x / Nx;                    % [m]
dy = dx;                        % [m]
dz = dx;                        % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
c0 = 1540;                      % [m/s]
rho0 = 1000;                    % [kg/m^3]
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;

% create the time array
t_end = (Nx * dx) * 2.2 / c0;   % [s]
kgrid.makeTime(c0, [], t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 2e6 / sc;     % [Hz] 
tone_burst_cycles = 3;

% create the input signal using toneBurst
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength ./ (c0 * rho0)) .* input_signal;
% figure;plot(input_signal);
% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% define the physical properties of the phased array transducer 
transducer.number_elements = 128 / sc;      % total number of transducer elements
transducer.element_width = 2;               % width of each element [grid points]
transducer.element_length = 2 / sc;        % length of each element [grid points]
transducer.element_spacing = 0;             % spacing (kerf  width) between the elements [grid points]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = c0;                    % sound speed [m/s]
transducer.focus_distance = 30e-3;              % focus distance [m]
transducer.elevation_focus_distance = 30e-3;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]
transducer.steering_angle_max = 40;             % maximum steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Hanning';
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = ones(transducer.number_elements, 1);

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

% making the sound and density maps into 3-dimensions
sound_speed_fullmap = zeros(Nx,Ny,Nz);
density_fullmap = zeros(Nx,Ny,Nz);
for i = 1:Nz
    sound_speed_fullmap(:,:,i) = sound_speed_map;
    density_fullmap(:,:,i) = density_map; 
end

medium.sound_speed = sound_speed_fullmap;
medium.density = density_fullmap;

% =========================================================================
% RUN THE SIMULATION
% =========================================================================
% range of steering angles to test
steering_angles = linspace(-40,40,101);
% preallocate the storage
number_scan_lines = length(steering_angles);
scan_lines = zeros(number_scan_lines, kgrid.Nt);
% set the input settings
input_args = {...
    'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};

% run the simulation if set to true, otherwise, load previous results
if  RUN_SIMULATION
    
    % loop through the range of angles to test
    for angle_index = 1:number_scan_lines
        
        % update the command line status
        disp('');
        disp(['Computing scan line ' num2str(angle_index) ' of ' num2str(number_scan_lines)]);
        
        % update the current steering angle
        transducer.steering_angle = steering_angles(angle_index);
        
        % run the simulation
        sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});
        rf_file = strcat('rf_',num2str(angle_index));
        save(rf_file,'sensor_data');
        % extract the scan line from the sensor data
        scan_lines(angle_index, :) = transducer.scan_line(sensor_data);
        
    end
    
    % save the scan lines to disk
    save example_us_phased_array_scan_lines scan_lines;
    
else
    
    % load the scan lines from disk, if RUN_SIMULATION = false, previous data
    % will load from this file
    load example_us_phased_array_scan_lines
    
end

% trim the delay offset from the scan line data
t0_offset = round(length(input_signal) / 2) + (transducer.appended_zeros - transducer.beamforming_delays_offset);
scan_lines = scan_lines(:, t0_offset:end);

% get the new length of the scan lines
Nt = length(scan_lines(1, :));

% =========================================================================
% PROCESS THE RESULTS
% =========================================================================

% -----------------------------
% Remove Input Signal
% -----------------------------

% create a window to set the first part of each scan line to zero to remove
% interference from the input signal
scan_line_win = getWin(Nt * 2, 'Tukey', 'Param', 0.05).';
scan_line_win = [zeros(1, t0_offset * 2), scan_line_win(1:end/2 - t0_offset * 2)];

% apply the window to each of the scan lines
scan_lines = bsxfun(@times, scan_line_win, scan_lines);

% -----------------------------
% Time Gain Compensation
% -----------------------------

% create radius variable
r = c0 * (1:Nt) * kgrid.dt / 2;    % [m]

% create time gain compensation function based on attenuation value,
% transmit frequency, and round trip distance
tgc_alpha = 0.2;       % [dB/(MHz cm)]
tgc = exp(2 * tgc_alpha * tone_burst_freq * 1e-6 * r * 100);

% apply the time gain compensation to each of the scan lines
scan_lines = bsxfun(@times, tgc, scan_lines);

% -----------------------------
% Frequency Filtering
% -----------------------------

% filter the scan lines using both the transmit frequency and the second
% harmonic
scan_lines_fund = gaussianFilter(scan_lines, 1/kgrid.dt, tone_burst_freq, 30, true);  
scan_lines_harm = gaussianFilter(scan_lines, 1/kgrid.dt, 2 * tone_burst_freq, 30, true);

% -----------------------------
% Envelope Detection
% -----------------------------

% envelope detection
scan_lines_fund = envelopeDetection(scan_lines_fund);
scan_lines_harm = envelopeDetection(scan_lines_harm);

% -----------------------------
% Log Compression
% -----------------------------

% normalised log compression
compression_ratio = 3;
scan_lines_fund = logCompression(scan_lines_fund, compression_ratio, true);
scan_lines_harm = logCompression(scan_lines_harm, compression_ratio, true);

% -----------------------------
% Scan Conversion
% -----------------------------

% set the desired size of the image
image_size = [Nx * dx, Ny * dy];

% convert the data from polar coordinates to Cartesian coordinates for
% display
b_mode_fund = scanConversion(scan_lines_fund, steering_angles, image_size, c0, kgrid.dt);
b_mode_harm = scanConversion(scan_lines_harm, steering_angles, image_size, c0, kgrid.dt);

% =========================================================================
% VISUALISATION OF ULTRASOUND IMAGE
% =========================================================================

% create the axis variables
x_axis = [0, Nx * dx * 1e3];    % [mm]
y_axis = [0, Ny * dy * 1e3];    % [mm]

R_val = linspace(0,0.06,length(scan_lines_fund));
Theta_val = sin(linspace(-40,40,101).*pi/180);
[R_Mat, Theta_Mat] = meshgrid(R_val, Theta_val); % Make grid
[x,y] = pol2cart(asin(Theta_Mat),R_Mat); % Convert to Cartesian coordinates

figure;h11=surf(y*1000,x*1000,scan_lines_harm); colormap(gray);set(h11,'LineStyle','none');
xt = get(gca, 'XTick');set(gca, 'FontName','Calibri','FontSize', 20, 'FontWeight', 'bold');
grid off;view(2);title('US Image','FontName','Calibri','FontSize',20,'FontWeight','bold');box on;colorbar;axis tight;

% interpolation code section
R_val_new = linspace(0,0.06,length(scan_lines_fund)*10); % interpolating by factor of 10
Theta_val_new = sin(linspace(-40,40,101*10).*pi/180);
[R_Mat_new, Theta_Mat_new] = meshgrid(R_val_new, Theta_val_new); % Make grid
[x_new,y_new] = pol2cart(asin(Theta_Mat_new),R_Mat_new); % Convert to Cartesian coordinates

interp_image = interp2(R_Mat,Theta_Mat,scan_lines_harm,R_Mat_new,Theta_Mat_new); % interpolated image

figure;h11=surf(y_new*1000,x_new*1000,20*log10(interp_image)); colormap(gray);set(h11,'LineStyle','none');
xt = get(gca, 'XTick');set(gca, 'FontName','Calibri','FontSize', 20, 'FontWeight', 'bold');
grid off;view(2);title('US Image','FontName','Calibri','FontSize',20,'FontWeight','bold');box on;colorbar;axis tight;

% filtering code section
filtered_image = filter2(fspecial('average',3),interp_image)/255; % filtered image

figure;h11=surf(y_new*1000,x_new*1000,20*log10(filtered_image)); colormap(gray);set(h11,'LineStyle','none');
xt = get(gca, 'XTick');set(gca, 'FontName','Calibri','FontSize', 20, 'FontWeight', 'bold');
grid off;view(2);title('US Image','FontName','Calibri','FontSize',20,'FontWeight','bold');box on;colorbar;axis tight;
