from math import sin, cos, asin, atan2, modf, pi, sqrt

# Sidereal Time

def current_GMST(now):
    # Method from TS Kelso's paper "Orbital Coordinate Systems Part II"
    # http://celestrak.com/columns/v02n02/
    # Also, this page: http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html
    # theta_g(tau) = theta_g(0h) + rotation_of_earth*seconds_of_day
    # theta_g(0h) = 24110.54841 + 8640184.812866 Tu + 0.093104 Tu ** 2 - 6.2 * 10e-6 Tu ** 3
    # Tu = du/36525
    # du = julian days since Jan 1, 2000 aka 2451545.0
    #
    current_JD  = current_julian_day(now)
    UT          = frac (current_JD + 0.5);
        # The purpose of UT is to flip it by 12h,
        # since Noon is the top of Julian Day
        # and Midnight is the top of UTC day
    current_JD  = current_JD - UT;
    du          = current_JD - 2451545.0 # Jan 1 2000
    Tu          = du/float(36525)
    GMST_0hour_seconds  = 24110.54841 + (Tu * (8640184.812866 + (Tu * (0.093104 - (Tu * 6.2E-6)))))
    GMST_0hour_seconds  = (GMST_0hour_seconds + (86400.0*1.00273790934*UT)) % 86400
    GMST_0hour_angle = (GMST_0hour_seconds/86400.0) * 2 * pi
    GMST = GMST_0hour_angle + (7.29211510e-5 * seconds_of_day(now))
    return GMST

def current_julian_day(now):
    current_jd = (julian_day_jan_1(now) + day_of_year(now))
    return current_jd

def julian_day_jan_1(now):
    # teh Julian Day of Jan 1st this year, first step in getting full julian day
    year = now.year - 1
    A = year / 100
    B = (2 - A) + trunc(A / 4)
    return (trunc(365.25 * year) + trunc(30.6001 * 14) + 1720994.5 + B)

def day_of_year(now):
    # Count all days in current year up to current day
    year = now.year
    month = now.month
    day = now.day
    days_of_months = (31,28,31,30,31,30,31,31,30,31,30,31)
    num_of_days = 0
    month_itor = 0;
    while ((month_itor+1) < month):
        num_of_days += days_of_months[month_itor]
        month_itor += 1
    num_of_days += day
    if (((year % 4) == 0) and
        (((year % 100) != 0) or ((year % 400) == 0)) and
        (month > 2)):
        num_of_days += 1;
    return num_of_days

def seconds_of_day(now):
    current_second = ((now.hour)*60*60) + (now.minute*60) + now.second
    return current_second

# Coordinate Frame Transforms
def spherical(coords, gmst):
    Xe  = coords[0]
    Ye  = coords[1]
    Ze  = coords[2]
    Re  = a   = 6378.137
    radius                  = sqrt  (Xe**2 + Ye**2 + Ze**2) - Re
    longitude               = (atan2(Ye, Xe)) - gmst
    latitude                = (atan2(Ze, sqrt ((Xe**2) + (Ye**2))))
    return [longitude, latitude, radius]

# Coordinate Frame Transforms
def spherical(coords):
    Xe  = coords[0]
    Ye  = coords[1]
    Ze  = coords[2]
    Re  = a   = 6378.137
    radius                  = sqrt  (Xe**2 + Ye**2 + Ze**2) - Re
    longitude               = (atan2(Ye, Xe))
    latitude                = (atan2(Ze, sqrt ((Xe**2) + (Ye**2))))
    return [longitude, latitude, radius]

def eci_to_ecf(eci_coords, gmst):
    # ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
    #
    # [X]     [C -S  0][X]
    # [Y]  =  [S  C  0][Y]
    # [Z]eci  [0  0  1][Z]ecf
    #
    #
    # Inverse:
    # [X]     [C  S  0][X]
    # [Y]  =  [-S C  0][Y]
    # [Z]ecf  [0  0  1][Z]eci
    #
    #

    X = (eci_coords[0] * cos(gmst))    + (eci_coords[1] * sin(gmst))
    Y = (eci_coords[0] * (-sin(gmst))) + (eci_coords[1] * cos(gmst))
    Z =  eci_coords[2]
    return [X, Y, Z]

def ecf_to_eci(ecf_coords, gmst):
    # ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
    #
    # [X]     [C -S  0][X]
    # [Y]  =  [S  C  0][Y]
    # [Z]eci  [0  0  1][Z]ecf
    #
    X = (ecf_coords[0] * cos(gmst))    - (ecf_coords[1] * sin(gmst))
    Y = (ecf_coords[0] * (sin(gmst)))  + (ecf_coords[1] * cos(gmst))
    Z =  ecf_coords[2]
    return [X, Y, Z]

def eci_to_geodetic(eci_coords, gmst, initial_lat):
    # http://www.celestrak.com/columns/v02n03/
    Xe  = eci_coords[0]
    Ye  = eci_coords[1]
    Ze  = eci_coords[2]
    a   = 6378.137
    b   = 6356.7523142
    R   = sqrt((Xe**2) + (Ye**2))
    f   = (a - b)/a
    e2  = ((2*f) - (f**2))
    longitude = atan2 (Ye, Xe) - gmst
    kmax = 20
    k = 0
    latitude = initial_lat
    while (k < kmax):
        C = 1 / sqrt( 1 - ( (e2)*(sin(latitude)**2) ) )
        latitude = atan2 (Ze + (a*C*e2*sin(latitude)), R)
        k += 1
    h = (R/cos(latitude)) - (a*C)
    return [longitude, latitude, h]

def ecf_to_geodetic(ecf_coords):
    """
        Conversion of Earth-Centered Earth-Fixed to Geodetic Latitude and Longitude

        (I'm using this method, Wu-Wang-Hu method)
        Yuanxin Wu; Ping Wang; Xiaoping Hu; , "Algorithm of Earth-centered Earth-fixed coordinates to geodetic coordinates," Aerospace and Electronic Systems, IEEE Transactions on , vol.39, no.4, pp. 1457- 1461, Oct. 2003
        doi: 10.1109/TAES.2003.1261144
        keywords: {Approximation algorithms;Earth;Ellipsoids;Equations;Geodesy;Global Positioning System;Iterative algorithms;Iterative methods;Kernel;Navigation; Newton-Raphson method; geodesy; navigation; ECEF; Earth-centered Earth-fixed coordinate; GPS navigation; NR method; Newton-Raphson method; geodetic coordinate; iterative approach;}
        URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1261144&isnumber=28184

        Zhu, J.; , "Conversion of Earth-centered Earth-fixed coordinates to geodetic coordinates," Aerospace and Electronic Systems, IEEE Transactions on , vol.30, no.3, pp.957-961, Jul 1994
        doi: 10.1109/7.303772
        keywords: {Application software;Approximation methods;Coordinate measuring machines;Earth;Ellipsoids;Equations;Geodesy;Navigation;Roundoff errors;Transforms;geodesy;geophysics computing;radionavigation;satellite relay systems;transforms;ECEF coordinates;Earth-centered Earth-fixed coordinates;NAVSTAR/GPS navigation;complexity;computer round-off error;geodetic coordinates;sensitivity;transformation;}
        URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=303772&isnumber=7475

        Olson, D.K.; , "Converting Earth-centered, Earth-fixed coordinates to geodetic coordinates," Aerospace and Electronic Systems, IEEE Transactions on , vol.32, no.1, pp.473-476, Jan. 1996
        doi: 10.1109/7.481290
        keywords: {Equations;Geometry;Government;Microcomputers;Testing;US Government;geodesy;geophysical signal processing;subroutines;-3000000 to 3000000 m;Earth-centered Earth-fixed coordinates;coordinates conversion;coversion routine;geodetic coordinates;improved algorithm;testbed program;}
        URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=481290&isnumber=10254
    """
    """
        P is the ECEF Position of the satellite
        P = (Xe Ye Ze) ECEF Coordinates

        R   length of the semimajor axis
        R   = 6378.137 km based on Model WGS-84
        r   length of the semiminor axis
        r   = 6356.7523142 km based on Model WGS-84
        f   flattening of the ellipsoid
        f   = (R - r)/R
        For the point P
        L   geodetic latitude,
        EW  geodetic longitude,
        h   altitude normal to reference ellipsoid.


        Latitude is easy:
        EW = atan2(Ye, Xe)

        Cross Section view of the Earth,
        The Vertical Plane containing the Satellite (P),
        The Center of the Earth (O),
        and the Ze Axis (cuz it's vertical, duh)

        This Cross Section is an ellipse that can be represented as:
        (x^2/R^2) + (z^2/r^2) = 1
        x and z are the points on the ellipse (surface of earth)

        The coordinates of P in the OXZ frame (Cross Section view) are
        (x0, z0)
        which are given by
        x0 = sqrt (Xe^2 + Ye^2)
        z0 = Ze

        We can write the ellipse as:
        x = R*cos(theta)
        z = r*sin(theta)
        -pi/2 < theta < pi/2
        x >= 0

        P1 is a point on the surface of the ellipse such that
        the normal of the ellipse at P1 will intersect P (satellite).
        Basically, is the point above Earth that the satellite is
        directly overhead.

        P2 is where the Ellipse intersects the Line going from P to O.
        This PO line is not normal to the ellipse.

        Coordinates for P1 can be written as:
        P1 = ( R cos (theta), r sin (theta) )

        The vector normal to the ellipse at P1 is:
        P1_n = ( cos (theta) / R, sin (theta) / r )

        The geodetic Longitude is the angle
        between the normal and the XY plane (equator).

        tan L = (sin(theta)/r) / (cos(theta)/R)
        tan L = (sin(theta)/cos(theta)*(R/r)))
        tan L = (tan(theta))*(R/r)

        Or:
        cot L = (cos(theta)/sin(theta)) * (r/R)
        cot L =     cot(theta)          * (0 + r/R)
        cot L =     cot(theta)          * (1 - R/R) + r/R)
        cot L =     cot(theta)          * (1 - (R-r)/R) ;  f = (R-r)/R
        cot L = (1 - f) cot(theta)

        The line normal to the ellipse at point P1, which
        passes through P, is specified as
        (x0-R*cos(theta)) / (z0-r*sin(theta)) = r*cos(theta) / R*sin(theta)


        With a second parameter t = tan(theta/2), (6) can be
        rewritten in terms of the variable t. After rearranging
        the resultant equation, the parameter t is observed to
        satisfy the condition f(t) = 0 where:
        f(t) = At^4 + (B + C)t^3 + (B - C)t + A
        A = r*z0
        B = 2*R*x0
        C = 2(R^2 - r^2)

        f(t) is a quartic equation whose roots can be approximated
        by the powerful Newton-Raphson (NR) method
        f(t)    =  At^4  +  (B + C)t^3 + (B - C)t + A
        f'(t)   = 4At^3  + 3(B + C)t^2 + (B - C)

        t_(k+1) = t_k - f(t_k)/f'(t_k)   k = 0,1,2...

        More specifically:
        t = t0 - f(t0)/f'(t0)

        We need a reliable t0. For this purpose, we can use the
        t associated with P2, since P2 is the Geocentric Lat, and
        is quite close anyway.

        Now to get t0.
        P1 and P2 in the OXZ frame
        (x1, z1) and (x2, z2), respectively.
        Using tan(theta) = 2t/(1 - t^2)
        (r/R)*tan(theta0) = z1/x1 ~= z2/x2 ~= z0/x0

        which is equivalent to:
        z0*t^2 + 2(1 - f)x0*t - z0 = 0
        We can ignore the other root, since th we have a good t0.
        t0 = (1 - f)(x0/z0) + sgn(z0) * sqrt( 1 + [(1 - f)(x0/z0)]^2 )

        We then use this with the Newton Raphson method as described above.
        One Iteration is enough to provide sufficiently accurate t.

        Using tan(theta) = 2t/(1 - t^2)
        We plug t into this equation:
        tan L = (tan(theta))*(R/r)
        L = atan ( tan(theta) * R/r)
        L = atan ((2t / (1-t^2)) * R/r )

        And there you have it, L = Geodetic Latitude

        Height is given by:
        z0 = z1 + h sin(L)
        Using sin theta = 2t=(1 + t^2) we obtain the
        geodetic height through the equation

        h = sgn(L)*(z0 - r*(2t/(1 + t^2))) * sqrt ( 1 + [(1-f)(1-t^2)/2t]^2 )
    """
    """
        Given coordinates Xe, Ye, Ze, ellipsoid parameters
        R, r, f, iterative accuracy delta, and maximum limit of
        iteration kmax:
        1)  Start.
        2)  Compute x0, z0 and coefficients: A, B, C.
        3)  If z0 = 0, let L = 0, h = x0 - R and go to Step 10;
            else go to Step 4.
        4)  Compute initial value t0 for iteration and initialize k = 0.
        5)  Set k = k + 1.
        6)  Compute tk as the function of t(k-1)
        7)  If |tk - t(k-1)| <= delta
            go to Step 9; else go to Step 8.
        8)  If k < kmax, go to Step 5;
            else the algorithm does NOT converge, go to Step 12.
        9)  Determine geodetic latitude L.
        10) Determine geodetic longitude.
        11) Determine geodetic height h.
        12) End.
    """
    # Step 1
    delta = 10e-9   # accuracy tolerance
    # Step 2
    Xe  = ecf_coords[0]
    Ye  = ecf_coords[1]
    Ze  = ecf_coords[2]
    # WSG 84
    R   = 6378.137
    r   = 6356.7523142
    f   = (R - r)/R
    one_minus_f = (1-f)
    x0  = sqrt((Xe**2) + (Ye**2))
    z0  = Ze
    A   = r*z0
    B   = 2*R*x0
    C   = 2*((R**2) - (r**2))
    h   = 0
    L   = 0.0
    # Step 9
    EW = atan2 (Ye, Xe)
    # Step 3
    if (abs(z0) < delta):
        h = x0 - R
    else:
        # Step 4
        t0   = (-one_minus_f*(x0/z0)) + ( sign(z0) * sqrt(1 + ((one_minus_f*(x0/z0))**2)) )
        kmax = 20
        k    = 0
        t    = t0
        while (k < kmax):
            # Step 5
            k           += 1
            # Step 6
            t0_to_4      = t0**4
            t0_to_3      = t0**3
            t0_to_2      = t0**2

            f_of_t       =   (A*t0_to_4)    +   ((B+C)*t0_to_3) + ((B-C)*t0) + A
            f_prime_of_t = (4*A*t0_to_3)    + (3*(B+C)*t0_to_2) +  (B-C)

            t            = t0 - (f_of_t/f_prime_of_t)

            # Step 7
            if (abs(t - t0) < delta):
                break
            t0 = t
            # Step 8 - Loop again


        # L = acot2 ( ((1 - f)*(1 - (t**2))) , (2*t) )
        # L = atan2((2*t), ((1 - f)*(1 - (t**2))))
        t_to_2      = t**2
        t_times_2   = t*2

        # Step 10
        L = ( atan2(-t_times_2 , (one_minus_f)*(1-t_to_2)) )

        # Step 11
        parta   = ( z0 - (r*t_times_2/(1+t_to_2)) )
        partb   = sqrt( 1  + (((one_minus_f*(1-t_to_2)/t_times_2))**2) )
        # WTF
        h       = (sign(L) * parta * partb)
        #h       = sgn(L) * (z0 - r*(2t/(1 + t^2))) * sqrt ( 1 + [(1-f)(1-t^2)/2t]^2 )
    return [EW, L, h]

def geodetic_to_ecf(geodetic_coords):
    longitude   = geodetic_coords[0]
    latitude    = geodetic_coords[1]
    height      = geodetic_coords[2]
    a           = 6378.137
    b           = 6356.7523142
    f           = (a - b)/a
    e2          = ((2*f) - (f**2))
    normal      = a / sqrt( 1 - (e2*(sin(latitude)**2)))

    X           = (normal + height) * cos (latitude) * cos (longitude)
    Y           = (normal + height) * cos (latitude) * sin (longitude)
    Z           = ((normal*(1-e2)) + height) * sin (latitude)
    return [X, Y, Z]

def ecf_to_topocentric(observer_coords, satellite_coords):
    # http://www.celestrak.com/columns/v02n02/
    # TS Kelso's method, except I'm using ECF frame
    # and he uses ECI
    longitude   = observer_coords[0]
    latitude    = observer_coords[1]
    height      = observer_coords[2]

    observer_ecf = geodetic_to_ecf (observer_coords)

    rx      = satellite_coords[0] - observer_ecf[0]
    ry      = satellite_coords[1] - observer_ecf[1]
    rz      = satellite_coords[2] - observer_ecf[2]

    top_s   = ( (sin(latitude)*cos(longitude)*rx) + \
                (sin(latitude)*sin(longitude)*ry) - \
                (cos(latitude)*rz))
    top_e   = ( -sin(longitude) * rx) + (cos(longitude) * ry)
    top_z   = ( (cos(latitude)*cos(longitude)*rx) + \
                (cos(latitude)*sin(longitude)*ry) + \
                (sin(latitude)*rz))

    return [top_s, top_e, top_z]

def topocentric_to_look_angles(topocentric):
    top_s = topocentric[0]
    top_e = topocentric[1]
    top_z = topocentric[2]
    range_sat    = sqrt((top_s**2) + (top_e**2) + (top_z**2))
    El      = asin (top_z/range_sat)
    Az      = atan2 (-top_e, top_s) + pi
    return [Az, El, range_sat]

# Doppler
def doppler(position_ecf, velocity_ecf, my_location_ecf, listen_frequency):
    # http://www.dtic.mil/cgi-bin/GetTRDoc?AD=ADA472033
    current_range = sqrt((position_ecf[0] - my_location_ecf[0])**2 +
                         (position_ecf[1] - my_location_ecf[1])**2 +
                         (position_ecf[2] - my_location_ecf[2])**2)

    next_pos    = [position_ecf[0] + velocity_ecf[0],
                   position_ecf[1] + velocity_ecf[1],
                   position_ecf[2] + velocity_ecf[2]]
    next_range = sqrt((next_pos[0] - my_location_ecf[0])**2 +
                      (next_pos[1] - my_location_ecf[1])**2 +
                      (next_pos[2] - my_location_ecf[2])**2)
    range_rate = (next_range - current_range)

    v = sqrt((velocity_ecf[0]**2) +
             (velocity_ecf[1]**2) +
             (velocity_ecf[2]**2))

    v *= sign(range_rate)
    c = 299792.458
    f = ((c / (c + v)) * listen_frequency)
    return f

# Utility
def sign(value):
    # http://stackoverflow.com/questions/1986152/why-python-doesnt-have-a-sign-function
    if (value < 0):  return (-1)
    if (value >= 0): return 1

def trunc(value):
    frac_part, int_part = modf (value)
    return int_part

def frac(value):
    frac_part, int_part = modf (value)
    return frac_part

def degrees_long(radians):
    degrees = (radians/pi*180) % (360)
    if degrees > 180:
        degrees = 360 - degrees
        degrees = str (degrees) + "E"
    else:
        degrees = str (degrees) + "W"
    return degrees

def degrees_lat(radians):
    if radians > pi/2 or radians < (-pi/2): return "Err"
    degrees = (radians/pi*180)
    if degrees < 0:
        degrees = str (degrees) + "S"
    else:
        degrees = str (degrees) + "N"
    return degrees

