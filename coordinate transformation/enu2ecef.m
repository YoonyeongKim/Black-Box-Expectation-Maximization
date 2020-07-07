%convertion of EastNorthUP coordinates to ECEF coordinates
function [x] = enu2ecef(R_lat,R_long,R_az,a,b,c)

x = (inv([ cos(R_az) sin(R_az) 0;
          -sin(R_az) cos(R_az) 0;
                   0        0  1]*... 
         [           -sin(R_long)             cos(R_long)          0;
          -cos(R_long)*sin(R_lat) -sin(R_lat)*sin(R_long) cos(R_lat);
           cos(R_long)*cos(R_lat)  cos(R_lat)*sin(R_long) sin(R_lat)]))*[a;b;c];
end      
