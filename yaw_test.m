clear all;
clc;
addpath('.\functions\');


%% test case: Initial Heading = 10Deg, angular = 10Deg/s, dt = 1s, Post Heading = 20Deg
heading = 10
heading = 10/180*pi; %
q = zeros(4, 1);
q = euler2q(heading, 0, 0);
Cbn = q2dcm(q);
Cnb = Cbn';
gyro = [0, 0, 10/180*pi];
Wpbb = gyro;

dq = zeros(4, 1);
dq(1) = -(Wpbb(1)*q(2) + Wpbb(2)*q(3) + Wpbb(3)*q(4))/2;
dq(2) = (Wpbb(1)*q(1) + Wpbb(3)*q(3) - Wpbb(2)*q(4))/2;
dq(3) = (Wpbb(2)*q(1) - Wpbb(3)*q(2) + Wpbb(1)*q(4))/2;
dq(4) = (Wpbb(3)*q(1) + Wpbb(2)*q(2) - Wpbb(1)*q(3))/2;

dt = 1;
q = q + dq*dt; % use half quaternion to update heading
q = q_norm(q);
Cbn = q2dcm(q);
[yaw, pitch, roll] = dcm2euler(Cbn);
heading = yaw*180/pi