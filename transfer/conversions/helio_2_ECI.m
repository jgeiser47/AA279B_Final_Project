function C = helio_2_ECI()
% Convert from heliocentric ecliptic reference frame to ECI (equatorial)

obliquity = deg2rad(23.45);
C = [1,               0,               0; ...
     0, +cos(obliquity), -sin(obliquity); ...
     0, +sin(obliquity), +cos(obliquity)];
end