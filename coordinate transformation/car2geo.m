%conversion of cartesian coordinates to geographic coordinates
function [lat,long] = car2geo(x,y,z)


thet = 180/pi*acos(z/sqrt(x^2+y^2+z^2)) ;
phi = 180/pi*acos(x/sqrt(x^2+y^2))    ;


if z >= 0
    lat = 90 - thet ;
    else
        lat = thet - 90 ;
end
  
if y < 0
    long = 360 - phi ;
    else
        long = phi ;


end
