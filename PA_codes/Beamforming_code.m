% OVERVIEW: DELAY AND SUM BEAMFORMING FOR PA IMAGING
%
% Input is "sensor_data" output from the back-solver code
% Output is the beamformed image

clc;clear all;close all;

n_wave = 1; % USER INPUT number of wavelengths
n_sources = 128;% USER INPUT
wv_vect = 800; % USER INPUT for Wavelengths needed
dt = 4e-8; % USER INPUT as per sampling rate
for i = 1:n_wave
    load(strcat('sensor_data',num2str(wv_vect(i)),'.mat'));     
    sensor_data_temp=sensor_data_temp1;
    %------------------------Remove the samples--------------------------------
    NUM_ZEROED_SAMPLES = 30;
    NUM_ZEROED_SAMPLES_END =0;
    sensor_data_temp(:,1:NUM_ZEROED_SAMPLES) = 0;
    sensor_data_temp(:,end-NUM_ZEROED_SAMPLES_END:end) = 0;    
    
    C=1500; % Speed of sound in medium [m/s]
    [row,col] = size(sensor_data_temp);
        
    RF = sensor_data_temp';
    
    Fs=1/dt; % Sampling Rate [Hz]
    Num_Zeros=30; % Number of initial samples to zero
    
    %Transducer Parameters
    Num_Channels=128; % Number of active channels 
    Element_Pitch=200e-6; % Transducer array element pitch [m] USER INPUT
    
    % Imaging Parameters
    Num_Samples=length(sensor_data_temp); % Number of samples recorded for each element for each beam
    Num_Beams=102; % Number of beams transmitted for each PA frame
    R_Init=0.00; % Minimal radius [m]
    R_Res=0.00015; % Radial resolution [m]
    Num_R=396; %435 Number of radial grid elements
    Sin_Theta_X_Init=-0.699999988; % Minimal sin(angle)
    Sin_Theta_X_Res=0.014; % Angular resolution (actually resolution of Sin(Theta_X) )
    Normalize=1; % 0 - Raw Image, 1 Log compression
    Dynamic_Range=70; % db for log compression
    
    Recon_PA=0; % 0 for PA
    PA_DELAY_OFFSET =0;
    % Generate Grid
    L_Transducer=Num_Channels*Element_Pitch; % Transducer array length
    x_Transducer=Element_Pitch/2+(-L_Transducer/2:Element_Pitch:L_Transducer/2-Element_Pitch); % Location of each array element
    
    R=R_Init+(0:R_Res:(R_Res*(Num_R-1))); % Radial axis
    Sin_Theta_X=Sin_Theta_X_Init+(0:Sin_Theta_X_Res:Sin_Theta_X_Res*(Num_Beams-2)); % Polar axis (again actually Sin(Theta_X) )
    
    [R_Mat, Sin_Theta_X_Mat] = meshgrid(R, Sin_Theta_X); % Make grid
    R_Mat=R_Mat'; Sin_Theta_X_Mat=Sin_Theta_X_Mat'; % Make matrices be in the currect orientation
    
    % Calculate Delay Matrix
    for ch=1:Num_Channels
        Dist_PA(:,:,ch)=Recon_PA*R_Mat+sqrt(x_Transducer(ch)^2+R_Mat.^2-2*x_Transducer(ch).*R_Mat.*Sin_Theta_X_Mat); % The recive+transmit distance
    end
    
    Delay_PA=round(Fs*Dist_PA/C+1); % Distance to delay in # samples
    BeamDelay_PA=zeros(Num_R,Num_Beams-1); % Photoacoustics uses only one beam
    
    % Filtering, Windowing & Apodization
    RF(1:Num_Zeros,:)=0; % Zero the first Num_Zeros Samples
    
    % Hilbert Transform
    [RF_Analytic, lower] = envelope(RF);
    RF_Analytic = abs(RF_Analytic);

    % Beamforming
    Image_PA=zeros(Num_R,Num_Beams-1);
    BeamNum_PA=Num_Beams; % Photoacoustics uses only one beam
  
    Image_PA_delayed = zeros(Num_R,Num_Channels,Num_Beams-1);
   
    for Ch=1:Num_Channels
        RF_Ch_PA=RF_Analytic(:,Ch); 
        index = Delay_PA(:,:,Ch)+BeamDelay_PA + PA_DELAY_OFFSET;
        Image_PA_delayed(:,Ch,:) = RF_Ch_PA(index);
    end

    apodiz_win = gausswin(Num_Channels*2,0.5);

    start_pt = floor(linspace(Num_Channels,1,Num_Beams-1));
    Image_PA11 = zeros(Num_R,Num_Beams-1);
    for apo = 1:Num_Beams-1
        index = start_pt(apo);
        apo_win = apodiz_win(index:index+Num_Channels-1);
        Image_PA11(:,apo) = Image_PA_delayed(:,:,apo)*apo_win;
    end
    [x,y] = pol2cart(asin(Sin_Theta_X_Mat),R_Mat); % Convert to Cartesian coordinates

    figure;h11=surf(y*1000,x*1000,double(Image_PA11)); colormap(hot);set(h11,'LineStyle','none');
    xt = get(gca, 'XTick');set(gca, 'FontName','Calibri','FontSize', 20, 'FontWeight', 'bold');    
    grid off;view(2);title('Linear PA Image','FontName','Calibri','FontSize',20,'FontWeight','bold');box on;colorbar;axis tight;
    
    figure;h11=surf(y*1000,x*1000,double(20*log10(Image_PA11))); colormap(hot);set(h11,'LineStyle','none');
    xt = get(gca, 'XTick');set(gca, 'FontName','Calibri','FontSize', 20, 'FontWeight', 'bold');    
    grid off;view(2);title('Log PA Image','FontName','Calibri','FontSize',20,'FontWeight','bold');box on;colorbar;axis tight;caxis([-10 0]);
end