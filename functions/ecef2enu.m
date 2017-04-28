function[enu] = ecef2enu(lla, ecef)

addpath('.\functions\');

sinlat = sin(lla(1));
coslat = cos(lla(1));
sinlon = sin(lla(2));
coslon = cos(lla(2));

enu(2) = -sinlat * coslon * ecef(1) - sinlat * sinlon * ecef(2) + coslat * ecef(3);
enu(1) = -sinlon * ecef(1) + coslon * ecef(2);
enu(3) = -(-coslat * coslon * ecef(1) - coslat * sinlon * ecef(2) - sinlat * ecef(3));
