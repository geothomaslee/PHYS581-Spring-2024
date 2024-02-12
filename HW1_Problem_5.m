%% 5a/b

clear all

alpha = 200; % Meters
fill_rate = 400000; % Cubic Meters Per Day

depth_per_day = fill_rate / (pi * alpha.^2); % Meters Per Day
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

%% 5c

v = 0.25; % Poisson's Ratio, Unitless
G = 10 * (1000^2); % Shear Modulus, Pascals (10 GPa)
Pz = pressure(24); % Pressure after 1 day of filling

% Displacement r < a
radius_steps = [0:5:alpha*5];
displacement = zeros(3,length(radius_steps));

for i = 1:length(radius_steps)
    radius = radius_steps(i); % current value of radius

    if i == 1
        displacement(1,i) = 0;
        displacement(2,i) = 0;
    elseif i > 1
        % rx and ry are the displacements, where r is the radius of a point
        % at position (x,y), and alpha is the radius of the load. X and Y
        % could be any (x,y), but we're only interested in abs(x) = abs(y)
        % which would be points along at a radius of alpha, which would 
        % have x and y coordinates of alpha / sqrt(2). Therefore, rx and ry
        % actually boil down to alpha^2 / sqrt(2)

        if radius <= alpha % Beneath Cylindrical Load, See Eq. 8.12
            scaling_factor = Pz * (1 - (2*v)) / (4 * G * alpha);

            displacement(1,i) = scaling_factor*(radius^2) * (1/sqrt(2)); % rx, but x in this case we only care about assuming radial movement from the center, so alpha
            displacement(2,i) = scaling_factor*(radius^2) * (1/sqrt(2)); % same as above for ry
    
            a = [0.5, -0.5];
            b = [1];
            z = (radius^2) / (alpha^2);
            hypergeom_func_comp = real(hypergeom(a, b, z));

            displacement(3,i) = scaling_factor * 4 * (alpha^2) * (1-v) * hypergeom_func_comp / (1-2*v);
        elseif radius > alpha % Away from Cylindrical Load, See Eq. 8.12
            scaling_factor = Pz * (alpha^2) * (1 - 2*v) / (4*G);

            displacement(1,i) = scaling_factor / (radius * sqrt(2));
            displacement(2,i) = scaling_factor / (radius * sqrt(2));

            a = [0.5, 0.5];
            b = [2];
            z = (alpha^2) / (radius^2);
            hypergeom_func_comp = real(hypergeom(a, b, z));

            displacement(3,i) = scaling_factor * 2 * (1-v) * hypergeom_func_comp / ((1-2*v) * radius);
        end
    end
end

r_alpha_ratio =  radius_steps/alpha; % After Figure 8.4

figure(3)
clf
hold on
plot(r_alpha_ratio,displacement(1,:),'Color','Blue','LineWidth',1)
plot(r_alpha_ratio,displacement(3,:),'Color','Red','LineWidth',1)
legend('Horizontal Displacement','Vertical Displacement','Location','Southeast')
xlabel('r / \alpha')
ylabel('Displacement')
title('Displacement from a Cylindrical Load of Radius \alpha')

figure(4)
clf
hold on
plot(radius_steps,displacement(1,:),'Color','Blue','LineWidth',1)
plot(radius_steps,displacement(3,:),'Color','Red','LineWidth',1)
legend('Horizontal Displacement','Vertical Displacement','Location','Southeast')
xlabel('Distance from Lake Center (meters)')
ylabel('Displacement')
title('Surface Displacements from Pu’u’Ō’ō Lava Lake (Meters)')
subtitle('After Filling At 400,000 Cubic Meters Per Day, for One Day')