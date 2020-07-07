function [ y ] = dyn_ballisticM_HV61_Measurement_YY( x , time, deltaTime, paramH)

global rx ry rz

%% Radar Position

px = x(1);
py = x(2);
pz = x(3);

R   = sqrt((px-rx)^2 + (py-ry)^2 + (pz-rz)^2);

Azm = atan2((py-ry),(px-rx));
Elv = atan2((pz-rz),sqrt((px-rx)^2+(py-ry)^2));

y = [R;
     Azm;
     Elv];
     
end