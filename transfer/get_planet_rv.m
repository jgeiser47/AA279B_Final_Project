function rv = get_planet_rv(ID, MJD)
    % Author: Joshua Geiser
    % Inputs: Planet ID number and Julian Date
    % Outputs: 6x1 array of heliocentric pos/vel [rx; ry; rz; vx; vy; vz]
    
    % Get orbital elements and rates at J2000 epoch
    [elems, rates] = AA279j2000_planetary_elements(ID);
    MJD_J2000 = cal_to_MJD(2000, 01, 01, 12, 0, 0);

    % Calculate JPL orbital elements at JD
    vals = zeros(6,1);
    for ii = 1:6
        vals(ii) = elems(ii) + rates(ii) * ((MJD - MJD_J2000) / 36525);
    end

    % Get standard orbital elements at JD
    a = au2km(vals(1));
    e = vals(2);
    i = wrapToPi(deg2rad(vals(3)));
    RAAN = wrapToPi(deg2rad(vals(4)));
    AOP = wrapToPi(deg2rad(vals(5)) - RAAN);
    M = wrapToPi(deg2rad(vals(6) - vals(5)));

    % Return helicoentric position and velocity
    rv = oe_to_rv(a, e, i, RAAN, AOP, M);
end