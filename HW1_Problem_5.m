%% 5a/b

clear all
clf

radius = 200; % Meters
fill_rate = 400000; % Cubic Meters Per Day

depth_per_day = fill_rate / (pi * radius.^2); % Meters Per Day
depth_per_hour = depth_per_day / 24; % Meters Per Hour
depth_initial = 0; % Meters

figure(1)
time = [1:1:24];
depth = (depth_per_hour * time) + depth_initial;

plot(time, depth)
title('Depth of Lava Lake')
subtitle('Given Initial Depth of 0')
xlabel('Hours Past Start of Day')
ylabel('Depth (m)')

% 5b
density = 2500; % kg/m^3
gravity = 9.81; % m/s
pressure = density * gravity * depth * -1; % Pascals
pressure_mpa = pressure/1000;

figure(2)
plot(time, pressure_mpa)
title('Pressure at Bottom of Lava Lake')
subtitle('Given Initial Depth of 0, -Pz is Downward')
xlabel('Hours Past Start of Day')
ylabel('Pressure (MPa)')