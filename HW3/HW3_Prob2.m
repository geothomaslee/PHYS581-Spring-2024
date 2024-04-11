%% Question 2)
clear all

V_steady = 0.1; % Initial velocity (m/s)
mu_steady = 0.6; % Initial frictional constant
a = 0.005;
b = 0.007;
Dc = 1 / 100; % Characteristic slip distance (m)
sigma = 100 * 10^6; % Normal stress (Pascals)
Fs_steady = mu_steady * sigma; % Initially, spring and friction are balanced, so initial friction force is equal to initial spring force
m_over_A=1*10^7; % kg/m^3
K_over_A = 1*10^9; % Pa/m
x_steady = 0; % Initial position of slider

psi_steady = Dc / V_steady; % Slip-rate

time_step = 0.01; % Time step in seconds
times = 0:time_step:10; % Times (seconds)

psi = zeros(1,length(times));
V = zeros(1,length(times));
Fs =  zeros(1,length(times));
x = zeros(1,length(times));
mu = zeros(1,length(times));
dpsi_dt = zeros(1,length(times));
V_load = 10*V_steady;

for t=1:length(times)
    if t == 1
        V(t) = V_steady; % at t=0, velocity suddenly jumps up to 10*V_steady

        psi(t) = V_steady / Dc; % Initially, the slider isn't moving
        
        dpsi_dt(t) = 1 - (V(t) * psi_steady / Dc);
        
        mu(t) = mu_steady + (a * log(V(t)/V_steady)) + (b*log(V_steady*psi(t) / Dc)); % Mu is calculated at this step
         
        x(t) = x_steady; % Even though Vload has gone up, at time = 0 it hasn't actually moved anything yet

        Fs(t) = Fs_steady; % The spring is unstretched after 0 time
    else
        psi(t) = psi(t-1) + dpsi_dt(t)*time_step; % Our new psi is calculated using the dpsi_dt we calculated in the last step
        x(t) = x_steady + (V_load*time_step) - (V(t-1)*time_step); % V_load stretches the string, but as the slider moves it will also change the stretch
        Fs(t) = K_over_A*x(t); % The spring has been stretched by some amount
        frictional_force = mu(t-1) * sigma; % Our current frictional force is estimated from the previous time stamp's coefficient of friction

        % If our spring force is stronger than the friction, our slider
        % will accelerate
        if Fs(t) > frictional_force
            V(t) = V(t-1) +((time_step/m_over_A)*(Fs(t) - frictional_force)); % Our new velocity is calculated from the difference in force
        else
            V(t) = V(t-1); % If not, the velocity just remains constant
        end
 
        mu(t) = mu_steady + (a * log(V(t)/V_steady)) + (b*log(V_steady*psi(t) / Dc)); % Now that the slider's velocity 
        % is changing, so should its coefficient of friction

    end
end

figure(1), clf
plot(times,mu)

