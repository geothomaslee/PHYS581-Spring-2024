%% Question 2)
clear all

V_steady = 1 / 100; % Initial velocity (m/s)
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

run_time = 30
dt = 0.01; % Time step in seconds
times = 0:dt:run_time; % Times (seconds)

psi = zeros(1,length(times));
V = zeros(1,length(times));
Fs =  zeros(1,length(times));
x = zeros(1,length(times));
mu = zeros(1,length(times));
dpsi_dt = zeros(1,length(times));
V_load = 10*V_steady;

for t=1:length(times)
    if t == 1
        V(t) = V_steady; % at t=0, the slider hasn't accelerated yet, it has the same initial velocity

        psi(t) = V_steady / Dc; % Initially, the slider isn't moving
        
        dpsi_dt(t) = 1 - (V(t) * psi(t) / Dc);
        
        mu(t) = mu_steady + (a * log(V(t)/V_steady)) + (b*log(V_steady*psi(t) / Dc)); % Mu is calculated at this step
         
        x(t) = x_steady; % Even though Vload has gone up, at time = 0 it hasn't actually moved anything yet

        Fs(t) = Fs_steady; % The spring is unstretched after 0 time
    else
        frictional_force = mu(t-1) * sigma; % Our current frictional force is estimated from the previous time stamp's coefficient of friction
        % If our spring force is stronger than the friction, our slider will accelerate
        if Fs(t-1) > frictional_force
            % Do nothing
        else
            frictional_force = Fs(t-1); % If not, the frictional force just matches the spring force (avoid backwards acceleration)
        end

        V(t) = V(t-1) + (dt/m_over_A)*(Fs(t-1) - frictional_force);

        x(t) = x(t-1) + (V_load - V(t-1))*dt;
        Fs(t) = Fs_steady + (K_over_A * x(t));

        psi(t) = psi(t-1) + (dpsi_dt(t-1) * dt);

        mu(t) = mu_steady + (a * log(V(t)/V_steady)) + (b*log(V_steady*psi(t) / Dc)); % Now that the slider's velocity 
        % is changing, so should its coefficient of friction

        dpsi_dt(t) = 1 - (V(t) * psi(t) / Dc);

    end
end

figure(1), clf

tt = (times*Dc) / V_steady

plot(tt,V-V_steady)


%% Testing functionn

dt = 0.01
times = 0:dt:4

steady_times = -2:1:0
steady_mu = mu_steady*ones(length(steady_times))

figure(2), clf
plot(steady_times,steady_mu,'k-')
hold on
K_over_As = [1e9:0.5e9:5e9]

for i = 1:length(K_over_As)
    K_over_A = K_over_As(i)
    friction = slider_function(times,dt,K_over_A)
    plot(times,friction)
    hold on
end





