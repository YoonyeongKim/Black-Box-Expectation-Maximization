function [ y ] = addRadarNoise(y_true,noiseLevel)

R = y_true(1);

Cr  = 0.000001    * 0.0001 * noiseLevel;   % Tuning parameter for the range measurement
Cth = 0.000000001 * 0.0001 * noiseLevel;   % Tuning parameter for the angle measurements

sigma_r  = Cr  * R^2;
sigma_th = Cth * R^2;

y = [y_true(1) + sigma_r*randn(1);
     y_true(2) + sigma_th*randn(1);
     y_true(3) + sigma_th*randn(1)];
     
end