function [B, V, W_inv, res] = cali4inv(Mag)

%% least square

% r = Y - X*beta
M = length(Mag); % measurement number
Y = zeros(M, 1);
X = zeros(M, 4);

for i = 1:M
    Bp = Mag(i, :);
    Y(i) = sum(Bp.*Bp);
    X(i, :) = [Bp, 1];
end

Beta = (X'*X)^-1*X'*Y;
Vx = 0.5*Beta(1);
Vy = 0.5*Beta(2);
Vz = 0.5*Beta(3);
V = [Vx, Vy, Vz];
W_inv = eye(3, 3);
B = sqrt(Beta(4) + Vx^2 + Vy^2 + Vz^2);

%% residual
E = (Y - X*Beta)'*(Y - X*Beta);
res = sqrt(E/M)/2/(B^2);

%% plot

% ellipsoid parameter
% (u - u0)'*A*(u - u0) = 1
A = eye(3, 3)/(B^2);
u0 = V;

% elipsode plot
if 1
figure;
% (u - u0)'*A*(u - u0) = 1
% A = R'* R
% sphere = R * (u - u0)
% u = inv(R) * sphere + u0
N_sphere = 50;
[x_sphere, y_sphere, z_sphere] = sphere(N_sphere);
% R = chol(A);
R = eye(3, 3)/B;
xyz_ellipsoid = inv(R)*[x_sphere(:)'; y_sphere(:)'; z_sphere(:)'] + u0' * ones(1, (N_sphere+1)*(N_sphere+1));
x_ellipsoid = reshape(xyz_ellipsoid(1,:), size(x_sphere));
y_ellipsoid = reshape(xyz_ellipsoid(2,:), size(y_sphere));
z_ellipsoid = reshape(xyz_ellipsoid(3,:), size(z_sphere));
surf(x_ellipsoid, y_ellipsoid, z_ellipsoid);
shading interp;
colormap(cool);
axis equal;
grid on;
xlabel('X');ylabel('Y');zlabel('Z');
hold on;

plot3(Mag(:,1), Mag(:, 2), Mag(:, 3), '*');
title(['4INV raw mag data: error = ', num2str(res*100), '%']);
end

% calibrated mag data
Mag_c = zeros(M, 3);
for i = 1:M
    Bp = Mag(i, :);
    Bc = Bp - V;
    Mag_c(i, :) = Bc;
end

if 0
[sphere_x, sphere_y, sphere_z]=sphere(50);
figure;
surf(B*sphere_x, B*sphere_y, B*sphere_z);
shading interp;
colormap(cool);
hold on;
plot3(Mag_c(:,1), Mag_c(:, 2), Mag_c(:, 3), '*');
title('4INV calibrated mag data');
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;
end