function[ecef] = lla2ecef(lla)

addpath('.\functions\');

RE = 6378137.0;
Rb = 6356752.3142;

EC2 = (RE * RE - Rb * Rb) / RE / RE;

sinlat = sin(lla(1));
coslat = cos(lla(1));
sinlon = sin(lla(2));
coslon = cos(lla(2));

N = RE / sqrt(1 - EC2 * sinlat * sinlat);
H = N + lla(3);

ecef(1) = H * coslat * coslon;
ecef(2) = H * coslat * sinlon;
ecef(3) = (N * (1 - EC2) + lla(3)) * sinlat;
