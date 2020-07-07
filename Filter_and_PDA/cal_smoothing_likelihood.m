function cal = cal_smoothing_likelihood(obs_enu, mu_x, variance_x, DY)

position_mu = mu_x(1:DY);
position_var = variance_x(1:DY,1:DY);

temp1 = (-1/2) * (obs_enu - position_mu)' * inv(position_var) * (obs_enu - position_mu);
temp2 = (-1/2) * log(det(position_var));
temp3 = (-1/2) * log((2*pi)^DY);

cal = temp1+temp2+temp3;

% cal = (1/sqrt((2*pi)^DY*det(position_var))) * exp((-1/2)*(obs_enu-position_mu)'*inv(position_var)*(obs_enu-position_mu));
     
end