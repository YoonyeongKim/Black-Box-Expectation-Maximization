%conversion of topographic coordinates to cartesian coordinates 
function [x, y, z] = top2car(az, el, latH, latD, lonH, lonD)
deg2rad = pi / 180;
if latH == 'N'
lat = latD  * deg2rad;
elseif latH == 'S'
lat = -latD * deg2rad;
end
if lonH == 'E'
lon = lonD * deg2rad;
elseif lonH == 'W'
lon = -lonD *deg2rad;
end
HA = sin(el);
EA = cos(el) * cos(az);
NA = cos(el) * sin(az);
%Rotation vector
T = [ -sin(lat)*cos(lon) -sin(lon) cos(lat)*cos(lon)
-sin(lat)*sin(lon) cos(lon) cos(lat)*sin(lon)
cos(lat) 0 sin(lat) ];
solVec = [0; 0; 0];
solVec = T * [EA; NA; HA];
x = solVec(1);
y = solVec(2);
z = solVec(3);
return