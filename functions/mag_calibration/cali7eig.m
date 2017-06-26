function [B, V, W_inv, res] = cali7eig(Mag)

%% least square

% r = X*beta
M = length(Mag); % measurement number
X = zeros(M, 7);

for i = 1:M
    Bp = Mag(i, :);
    X(i, :) = [Bp.^2, Bp, 1];
end

[eig_V, eig_D] = eig(X'*X);
E = eig_D(1, 1);
Beta = eig_V(:, 1);

A = diag(Beta(1:3));
if det(A) < 0
    A = -A;
    Beta = -Beta;
end
Axx = A(1, 1);
Ayy = A(2, 2);
Azz = A(3, 3);

Vx = -Beta(4)/2/Beta(1);
Vy = -Beta(5)/2/Beta(2);
Vz = -Beta(6)/2/Beta(3);
V = [Vx, Vy, Vz];

B = sqrt(abs(Axx*Vx^2 + Ayy*Vy^2 + Azz*Vz^2 - Beta(7)));

%% residual
res = sqrt(E/M)/2/(B^2);

%% normalize the ellipsoid matrix A to unit determinant
A_det = det(A);
A = A / A_det^(1/3);
% (Bp - V)'*A*(Bp - V) = B^2
B = B / A_det^(1/6);

%% plot

% ellipsoid parameter
% (Bp - V)'*A*(Bp - V) = (Bp - V)'*(W^-1)'*(W^-1)*(Bp - V) = B^2
W_inv = sqrt(A);
u0 = V;

% elipsode plot
if 0
figure;
% (u - u0)'*A*(u - u0) = 1
% A = R'* R
% sphere = R * (u - u0)
% u = inv(R) * sphere + u0
N_sphere = 50;
[x_sphere, y_sphere, z_sphere] = sphere(N_sphere);
R = W_inv/B;
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
alpha(0.3);
hold on;

plot3(Mag(:,1), Mag(:, 2), Mag(:, 3), '*');
title(['7EIG raw mag data: error = ', num2str(res*100), '%']);
end

% calibrated mag data
Mag_c = zeros(M, 3);
for i = 1:M
    Bp = Mag(i, :);
    Bc = W_inv*(Bp - V)';
    Mag_c(i, :) = Bc;
end

if 0
[sphere_x, sphere_y, sphere_z]=sphere(50);
figure;
surf(B*sphere_x, B*sphere_y, B*sphere_z);
shading interp;
colormap(cool);
alpha(0.3);
hold on;
plot3(Mag_c(:,1), Mag_c(:, 2), Mag_c(:, 3), '*');
title('7EIG calibrated mag data');
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;
end
