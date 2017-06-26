function[heading] = heading_between_2points(x1_lla, x2_lla)

addpath('.\functions\');

x1_latitude = x1_lla(1);
x1_longitude = x1_lla(2);
x2_latitude = x2_lla(1);
x2_longitude = x2_lla(2);

RE = 6378137.0;
esqu = 0.00669437999013;
ave_latitude = (x1_latitude + x2_latitude)/2;
RM = RE * (1 - esqu) / ((1 - esqu * sin(ave_latitude) * sin(ave_latitude))^1.5);
RN = RE / (sqrt(1 - esqu * sin(ave_latitude) * sin(ave_latitude)));

delta_latitude = x2_latitude - x1_latitude;
delta_longitude = x2_longitude - x1_longitude;

delta_N = delta_latitude * RM;
delta_E = delta_longitude * RN * cos(ave_latitude);

heading = atan2(delta_E, delta_N);