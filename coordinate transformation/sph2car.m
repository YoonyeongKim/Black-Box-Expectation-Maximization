%conversion of spherical coordinates to cartesian coordinates
function [x, y, z] = sph2car(thet, phi, R)

    x = sin(thet) * cos(phi) * R    ;
    y = sin(thet) * sin(phi) * R    ;   
    z = cos(thet) * R               ;

return