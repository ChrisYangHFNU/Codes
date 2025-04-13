% ================== Plot Vertical Coverage Pattern Using Default Parameters ================== 
% Set the frequency to 100MHz, the antenna height to 10m, and the free-space range to 200km,
% The antenna pattern, surface roughness, antenna tilt angle, and field
% polarization assume their default values as specified in the
% AntennaPattern, SurfaceRoughness, TiltAngle, and Polarization properties.
% Obtain an array of vertical coverage pattern values and angles
% 
% freq: Radar frequency(real-valued scalar less than 10 GHz)
% rfs: Free-space range(positive scalar|positive vector)
% anht: Radar antenna height(real-valued scalar)
clear all;
close all;

freq = 100e6;
ant_height = 10;
rng_fs = [100 200]; % unit:km
[vcp,vcpangles] = radarvcd(freq,rng_fs,ant_height);

% To see the vertical coverage pattern, omit the output arguments.
radarvcd(freq,rng_fs,ant_height);

% ================== Vertical Coverage Pattern with Specified Antenna Pattern ==================
% set the frequency to 100MHz, the antenna height to 10m, and the
% free-space range to 200km. The antenna pattern is a sinc function with
% 45° half-power width. The surface height standard deviation is set to
% 1/2sqrt(2)m. The antenna tilt angle is set to 0°, and the field
% polarization is horizontal
pat_angles = linspace(-90,90,361)';
freq = 100e6;

ntn = phased.SincAntennaElement('Beamwidth',45);
pat = ntn(freq,pat_angles');

ant_height = 10;
rng_fs = 200;
tilt_ang = 0;
[vcp,vcpangles] = radarvcd(freq,rng_fs,ant_height,...
    'RangeUnit','km','HeightUnit','m',...
    'AntennaPattern',pat,...
    'PatternAngles',pat_angles,...
    'TiltAngle',tilt_ang,'SurfaceHeightStandardDeviation',1/(2*sqrt(2)));

% Call radarvcd with no output arguments to display the vertical coverage
% pattern.
radarvcd(freq,rng_fs,ant_height,...
    'RangeUnit','km','HeightUnit','m',...
    'AntennaPattern',pat,...
    'PatternAngles',pat_angles,...
    'TiltAngle',tilt_ang,'SurfaceHeightStandardDeviation',1/(2*sqrt(2)));

% Alternatively, use the radarvcd output arguments and the blakechart
% function to display the vertical coverage pattern to a maximum range of
% 400 km and a maximum height to 50 km. Customize the Blake chart by
% changing the color.
blakechart(vcp,vcpangles,400,50,...
    'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor',[0.8500 0.3250 0.0980])

