function [obs_enu] = dyn_ballisticM_HV61_Measurement_YY_inv(obs_sph)

global rx ry rz

%% Radar Position

[tempx tempy tempz] = sph2cart(obs_sph(2), obs_sph(3), obs_sph(1));

enu_x = tempx + rx;
enu_y = tempy + ry;
enu_z = tempz + rz;

obs_enu = [enu_x; enu_y; enu_z];

     
end