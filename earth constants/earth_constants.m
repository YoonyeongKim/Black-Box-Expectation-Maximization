    

    % Earth Constants
    constant.posEarth = [0; 0; 0]                                  ; %Earth's center position
    constant.Re = 6371e3                                           ; %Earth radius (m)
    constant.Me = 5.9742e24                                        ; %Earth mass (kg)
    constant.Gc = 6.673e-11                                        ; %Gravitational constant (m^3 kg^-1 s^-2)
    constant.g0 = constant.Gc * constant.Me / (constant.Re ^ 2)    ; %Gravitational Acceleration (sea level)
    constant.c = 299792458                                         ; %Speed of light (m/s)
    constant.omega = [0 0 7.2921159e-5]'                           ; %angular speed (rad/s)