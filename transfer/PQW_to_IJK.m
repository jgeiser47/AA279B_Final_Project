function R_PQW_to_IJK = PQW_to_IJK(a, e, i, RAAN, AOP, MA)
    % Author: Joshua Geiser
    % Inputs: Keplerian elements (km, rad, rad, rad, rad, rad)
    % Output: DCM from perifocal to inertial coordinates

    % Circular, equatorial orbit
    if e == 0 && i == 0
        R_PQW_to_IJK = eye(3); 
        
    % Circular, non-equatorial orbit
    elseif e == 0
        R_PQW_to_IJK = [cos(RAAN), -sin(RAAN)*cos(i), +sin(RAAN)*sin(i);
                        sin(RAAN), +cos(RAAN)*cos(i), -cos(RAAN)*sin(i);
                        0        , +sin(i)          , +cos(i)];
        
    % Non-circular, equatorial orbit
    elseif i == 0
        R_PQW_to_IJK = [cos(AOP), -sin(AOP), 0;
                        sin(AOP), +cos(AOP), 0;
                        0,         0       , 0];
        
    % Non-circular, non-equatorial orbit (general case)
    else
        R_PQW_to_IJK = [cos(RAAN)*cos(AOP) - sin(RAAN)*cos(i)*sin(AOP), ...
                        -cos(RAAN)*sin(AOP) - sin(RAAN)*cos(i)*cos(AOP), ...
                        sin(RAAN)*sin(i); ...
                        sin(RAAN)*cos(AOP) + cos(RAAN)*cos(i)*sin(AOP), ...
                        -sin(RAAN)*sin(AOP) + cos(RAAN)*cos(i)*cos(AOP), ...
                        -cos(RAAN)*sin(i); ...
                        sin(i)*sin(AOP), ...
                        sin(i)*cos(AOP), ...
                        cos(i)];
    end
end