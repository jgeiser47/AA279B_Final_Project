function rv = oe_to_rv(a, e, i, RAAN, AOP, M)
    % Author: Joshua Geiser
    % Inputs: Keplerian elements (km, --, rad, rad, rad, rad)
    % Outputs: 3x2 matrix containing pos/vel vectors in inertial frame

    % Assumes heliocentric orbit
    mu = 1.3271244004193938e11;
    tol = 1e-10;

    % Calculate a couple necessary parameters
    E = MA_to_EA(M, e, tol);
    n = sqrt(mu / a^3);
    
    % Rotation matrix from perifocal to inertial (ECI) coordinates
    R_PQW_to_IJK = PQW_to_IJK(a, e, i, RAAN, AOP, M);
    
    % Position vector
    r_PQW = [a*(cos(E)-e); a*sqrt(1-e^2)*sin(E); 0];
    r_IJK = R_PQW_to_IJK * r_PQW;
    
    % Velocity vector
    v_PQW = ((a*n)/(1-e*cos(E))) .* [-sin(E); sqrt(1-e^2)*cos(E); 0];
    v_IJK = R_PQW_to_IJK * v_PQW; 
    
    rv = [r_IJK; v_IJK];
end