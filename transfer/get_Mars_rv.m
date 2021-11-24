function rv = get_Mars_rv(MJD)

model = 'ephemeris';

if strcmp(model, 'ephemeris')
    rv = get_planet_rv(4, MJD);
elseif strcmp(model, 'circular')
    
    % Assumes heliocentric orbit
    mu_Sun = 1.3271244004193938e11;
    
    % Ephemeris values for Mars on 01/28/2034 00:00:00 (arrival epoch)
    a = 227953237.299835; 
    e = 0.0934209570652977;
    i = 0.0322348483513284;
    RAAN = 0.863237189337131;
    AOP = -1.27848943960252;
    M = 1.06623161056081;

    % New values for a circular orbit
    a = 227953016; 
    n = sqrt(mu_Sun/(a^3));
    e = 0;
    
    % Reference values for mean motion and MJD at specific epoch
    M0 = M; 
    T0 = cal_to_MJD(2034, 1, 28, 0, 0, 0);
    
    M_curr = M0 + n * ((MJD - T0)*86400);
    rv = oe_to_rv(a, e, i, RAAN, AOP, M_curr);
end

err_str = sprintf("Mars propagation model must be set to 'ephemeris' or 'circular'\nGiven value of '%s' is invalid.", model);
assert((strcmp(model, 'ephemeris') || strcmp(model, 'circular')), ...
       err_str);
end