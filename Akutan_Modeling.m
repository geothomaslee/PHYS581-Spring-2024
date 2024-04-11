%% Mogi Model
step = 1; % Distance step
v = 0.25; % Poisson's Ratio
d = 10; % depth
mu = 1*(10^9); % Shear modulus
a = 10; % radius of sphere
V = 0.1;

p = (V * mu) / (pi * a^3);

[X,Y] = meshgrid(-10:step:10,-10:step:10);

Uz = zeros(length(X),length(Y));
Ux = zeros(length(X),length(Y));
Uy = zeros(length(X),length(Y));
 

for i = 1:length(X(:,1))
    for j = 1:length(Y(1,:))
        x = X(i,j);
        y = Y(i,j);
        r = sqrt(x^2 + y^2);
        Uz(i,j) = ((1-v)*p*(a^3) * d) / (mu * (r^2 + d^2)^(3/2)) *1000 * 100 * 10;

        Ur = ((1-v)*p*(a^3) * r) / (mu * (r^2 + d^2)^(3/2));
        if x < 0
            Ux(i,j) = Ur / sqrt(2) * -1 * 1000 * 100 * 10;
        else
            Ux(i,j) = Ur / sqrt(2) * 1000 * 100 * 10;
        end

        if y < 0
            Uy(i,j) = Ur / sqrt(2) * -1 * 1000 * 100 * 10;
        else  
            Uy(i,j) = Ur / sqrt(2) * 1000 * 100 * 10;
        end

    end
end

figure(3)
mesh(X+Ux,Y+Uy,Uz)
zlabel('Vertical Deformation (mm)')
xlabel('Kilometers')
ylabel('Kilometers')




