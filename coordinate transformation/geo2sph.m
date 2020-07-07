%conversion of geographic coordinates to spherical coordinates
function [thet, phi] = geo2sph(latH, latD , lonH, lonD)

deg2rad = pi / 180              ;
latDegree = latD                ;
latRad = latDegree * deg2rad    ;
lonDegree = lonD                ;
lonRad = lonDegree * deg2rad    ;

if latH == 'N'
    thet = pi / 2 - latRad      ;
    elseif latH == 'S'
        thet = pi / 2 + latRad  ;
end

if lonH == 'E'
    phi = lonRad                ;
    elseif lonH == 'W'
        phi = 2 * pi - lonRad   ;
end

return
