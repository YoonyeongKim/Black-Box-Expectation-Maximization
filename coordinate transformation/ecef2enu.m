%convertion of EarthCenteredEarthFixed coordinates to EastNorthUP coordinates
function [a,b,c] = ecef2enu(lat,long,az,x,y,z)
   
    % long & lat : radar position 
    a1 =          -sin(long)*x +          cos(long)*y                  ;
    b1 = -cos(long)*sin(lat)*x - sin(lat)*sin(long)*y + cos(lat)*z     ;
    c1 =  cos(long)*cos(lat)*x + cos(lat)*sin(long)*y + sin(lat)*z     ;
    
    % az : ENU(x=0 equals RadarAzmuthAngle=0)->specific AzAngle
    
    a =  cos(az)*a1 + sin(az)*b1    ; 
    b = -sin(az)*a1 + cos(az)*b1    ;
    c = c1                          ;
    
return



