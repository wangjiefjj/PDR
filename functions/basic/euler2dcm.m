function [cbn] = euler2dcm (yaw, pitch, roll)

cbn = zeros(3,3);

cbn(1, 1) = cos(pitch)*cos(yaw);
cbn(1, 2) = sin(roll)*sin(pitch)*cos(yaw) - cos(roll)*sin(yaw);
cbn(1, 3) = cos(roll)*sin(pitch)*cos(yaw) + sin(roll)*sin(yaw);

cbn(2, 1) = cos(pitch)*sin(yaw);
cbn(2, 2) = sin(roll)*sin(pitch)*sin(yaw) + cos(roll)*cos(yaw);
cbn(2, 3) = cos(roll)*sin(pitch)*sin(yaw) - sin(roll)*cos(yaw);

cbn(3, 1) = -sin(pitch);
cbn(3, 2) = sin(roll)*cos(pitch);
cbn(3, 3) = cos(roll)*cos(pitch);



