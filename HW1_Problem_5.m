%% 5a/b

clear all

alpha = 200; % Radius of the Load, Meters
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
pressure_mpa = pressure/1000000;

figure(2)
plot(time, pressure_mpa)
title('Pressure at Bottom of Lava Lake')
subtitle('Given Initial Depth of 0, -Pz is Downward')
xlabel('Hours Past Start of Day')
ylabel('Pressure (MPa)')

%% 5c
% Plotting the static deformation after a day of filling

v = 0.25; % Poisson's Ratio, Unitless
G = 10 * (10^9); % Shear Modulus, Pascals (10 GPa)
Pz = pressure(24); % Pressure after 1 day of filling

radius_steps = [0:5:alpha*5]; % Radius to step outwards
displacement = zeros(3,length(radius_steps));

num_iterations = length(radius_steps);

for i = 1:num_iterations
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

            displacement(3,i) = (scaling_factor * 4 * (alpha^2) * (1-v) * hypergeom_func_comp) / (1-2*v);
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
displacement_mm = displacement * 10 * 100;
plot(r_alpha_ratio,displacement_mm(1,:),'Color','Blue','LineWidth',1)
plot(r_alpha_ratio,displacement_mm(3,:),'Color','Red','LineWidth',1)
legend('Horizontal Displacement','Vertical Displacement','Location','Southeast')
xlabel('r / \alpha')
ylabel('Displacement (mm)')
title('Displacement from a Cylindrical Load of Radius \alpha')

figure(4)
clf
hold on
plot(radius_steps,displacement(1,:),'Color','Blue','LineWidth',1)
plot(radius_steps,displacement(3,:),'Color','Red','LineWidth',1)
legend('Horizontal Displacement','Vertical Displacement','Location','Southeast')
xlabel('Distance from Lake Center (Meters)')
ylabel('Displacement (mm)')
title('Surface Displacements from Pu’u’Ō’ō Lava Lake (Meters)')
subtitle('After Filling At 400,000 Cubic Meters Per Day, for One Day')

%% 5d
% Plotting the deformation over time

% Creating the model grid
X_range = -1000:50:1000;
Y_range = -1000:50:1000;
[X,Y] = meshgrid(X_range,Y_range);

days = 20; % Number of days to model
decay = 150; % A decay constant for the sinusoidal filling function
frequency = 0.1; % Frequency of the filling
time_steps = [0:1:days*24]; % Time steps, in hours
initial_depth = 0; % Initial depth of the lava lake
depth = zeros(1,length(time_steps)) + initial_depth; % Depth function, empty at first
inflow_base_rate = 400000; % Base inflow rate, m^3/s
inflow = 400000 * sin(frequency*time_steps).*exp(time_steps / -decay); % Inflow function
depth_change = inflow / (pi * alpha.^2); % Meters Per Hour

% Adding the depth change at each step to the previous value to find the
% current depth
for i = 2:length(depth)
    depth(i) = depth(i-1) + depth_change(i);
end

figure(7)
plot(time_steps,depth)
xlabel('Time (hours)')
ylabel('Depth (m)')
title('Model Sinusoidal Filling and Draining Lava Lake')

% Pressure, and max pressure
pressure = density * gravity * depth * -1;
max_pressure = max(abs(pressure));

% Calculating the maximum displacement using the max pressure at the center
a = [0.5, -0.5];
b = [1];
scaling_factor = max_pressure * (1 - (2*v)) / (4 * G * alpha);
z = 0;
max_displacement = scaling_factor * 4 * (alpha^2) * (1-v) * real(hypergeom(a, b, z)) / (1-2*v);

%% Break for calculation ease

% Creating a grid of the radius for every point, for the horizontal
% deformations
Radius = sqrt(X.^2 + Y.^2);
[rows, columns] = size(Radius);

% Creating a deformation grid temporarily just filled with zeros
Vertical_Deformation = zeros(rows,columns);

% Loops through every index (i,j) in the grid matrix and calculates a
% vertical deformation value

num_iterations = length(pressure); % Number of iterations, for testing. Minimum 2

for p = 2:num_iterations
    Vertical_Deformation = zeros(rows,columns);
    X_Deformation = zeros(rows, columns);
    Y_Deformation = zeros(rows, columns);
    
    Pz = pressure(p);

    for i = 1:rows
        for j = 1:columns
            if Radius(i,j) <= alpha % Near deformation
                a = [0.5, -0.5];
                b = [1];
                scaling_factor = Pz * (1 - (2*v)) / (4 * G * alpha);
                z = (Radius(i,j)^2) / (alpha^2);
                Vertical_Deformation(i,j) = scaling_factor * 4 * (alpha^2) * (1-v) * real(hypergeom(a, b, z)) / (1-2*v);
                X_Deformation(i,j) = -scaling_factor * Radius(i,j) * X(i,j); % Deformation in the X direction
                Y_Deformation(i,j) = -scaling_factor * Radius(i,j) * Y(i,j); % Deformation in the Y direction
            elseif Radius(i,j) > alpha % Distant deformation
                a = [0.5, 0.5];
                b = [2];
                scaling_factor = Pz * (alpha^2) * (1 - 2*v) / (4*G);
                z = (alpha^2) / (Radius(i,j)^2);
                Vertical_Deformation(i,j) = scaling_factor * 2 * (1-v) * real(hypergeom(a, b, z)) / ((1-2*v) * Radius(i,j));
                X_Deformation(i,j) = -scaling_factor * (X(i,j) / (Radius(i,j)^2)); 
                Y_Deformation(i,j) = -scaling_factor * (Y(i,j) / (Radius(i,j)^2));
            end
            % Replaces infinity values with 0
            if X_Deformation(i,j) > 1
                X_Deformation(i,j) = 0
            elseif Y_Deformation(i,j) > 1
                Y_Deformation(i,j) = 0
            end
        end
    end

    figure(5)
    clf(5)
    subplot(2,2,1)
    ylim auto, xlim auto
    surf(X,Y,Vertical_Deformation)
    colormap winter;
    title('Vertical Displacement')
    txt = ['\sigma_z = ' num2str(Pz/1000000) ' MPa'];
    subtitle(txt)
    zlabel('Displacement (m)')
    max_z = (-1.1 * max_displacement);
    zlim([max_z,0])

    subplot(2,2,2)
    max_depth = max(depth);
    height = 2 * max_depth;
    [X_cyl,Y_cyl,Z_cyl] = cylinder(alpha); % Cylinder
    [X_cyl_lava,Y_cyl_lava,Z_cyl_lava] = cylinder(alpha*1.002);
    Z_cyl_crater = Z_cyl * height;
    Z_cyl_lava = Z_cyl * depth(p);
    surf(X_cyl,Y_cyl,Z_cyl_crater,'FaceAlpha',0.2,'EdgeColor','none','FaceColor','Blue');
    hold on
    theta =0:0.01*pi:2*pi;
    x_lava_base=alpha*cos(theta); y_lava_base=alpha*sin(theta); z_lava_base=zeros(size(x_lava_base))+depth(p);
    % simple patch for base
    p_base =patch(x_lava_base, y_lava_base, z_lava_base,'Red','FaceAlpha',0.8,'EdgeColor',[1,0.02,0.02]);
    surf(X_cyl_lava,Y_cyl_lava,Z_cyl_lava,'FaceColor','Red','FaceAlpha',0.8,'EdgeColor',[1,0.02,0.02]);
    xlim([-1.1*alpha,1.1*alpha])
    ylim([-1.1*alpha,1.1*alpha])
    zlabel('Depth (m)')
    title(['Depth of Lava Lake: ' num2str(depth(p)) 'm'])
    subtitle(['Maximum Depth: ' num2str(max_depth) 'm'])
    
    subplot(2,2,3)
    quiver_scale = 100000;
    Scaled_X_Deformation = X_Deformation * quiver_scale;
    Scaled_Y_Deformation = Y_Deformation * quiver_scale;
    quiver(X,Y,Scaled_X_Deformation,Scaled_Y_Deformation,'off')
    xlim([-500,500])
    ylim([-500,500])
    xlabel('X (m)')
    ylabel('Y (m)')
    title('Horizontal Deformation')
    subtitle('Scaled by 10^5')

    subplot(2,2,4)
    deformation_scale = 2000 / max_displacement;
    Full_X_Def = X + (X_Deformation * deformation_scale);
    Full_Y_Def = Y + (Y_Deformation * deformation_scale);
    surf(Full_X_Def,Full_Y_Def,Vertical_Deformation)
    zlim([max_z,0])

    title('Full Deformation')
    subtitle('X/Y Displacements Exaggerated by 10^6')
    zlabel('Displacement (m)')

    % Supertitle for all the subplots
    R = mod(p, 24);
    Q = floor(p./24) + 1;
    txt = ['Deformation from Lava in Pu u Ōnō Lava Lake, Day ' num2str(Q) ' Hour ' num2str(R)];
    sgtitle(txt)
    
    filename = ['Vertical_Displacement_Animation.gif'];
    drawnow 
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    if p == 2; 
        imwrite(imind,cm,filename,'gif', 'DelayTime',0.02,'Loopcount',inf); 
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end
end

