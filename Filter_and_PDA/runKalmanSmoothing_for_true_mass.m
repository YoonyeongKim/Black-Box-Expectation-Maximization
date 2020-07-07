function smoothing_likelihood = runKalmanSmoothing_for_true_mass(mass0, delta_mass)

    global numMissile
    global startTime endTime timeLength deltaTime
    global startErrorLevel
    global DX DY beta_possible theta_possible y_obs ballistic_data
	global measure_cov movement_cov_initial movement_cov_final
    global alpha beta kappa gamma observationSize itrMissile
    
	%% Dynamics and Ballistic Data
    % Define dynamics
	F = 'dyn_ballisticM_HV61_thrust';
	F_inv = 'dyn_ballisticM_HV61_thrust_inv';
	H = 'dyn_ballisticM_HV61_Measurement_YY';
	H_inv = 'dyn_ballisticM_HV61_Measurement_YY_inv';

	% Define Mass of each time
	mass_initial = mass0;
	mass_final = mass_initial - timeLength*(deltaTime*10)*delta_mass;

	%% Forward filtering
    % Define initial state
	initial_x = ballistic_data.x(:,startTime);
	x_start = initial_x + initial_x*randn(1)*startErrorLevel;

	% Each filter represents UKF with candidate beta & candidate theta
    filters = cell(1);
	filterIMM = cell(numMissile,1);
    
	for b_idx = 1:length(beta_possible)
		for t_idx = 1:length(theta_possible)
			filters{3*(b_idx-1)+t_idx} = all_UKF_forward(F,[beta_possible(b_idx), theta_possible(t_idx)],H,[DY],x_start,movement_cov_initial,observationSize,alpha,beta,kappa,gamma,mass_initial);
		end
    end

	% Generate IMM integrating all the UKFs (filters)
	filterIMM{itrMissile} = all_IMM_forward(filters, x_start, movement_cov_initial, beta_possible, theta_possible, startTime, endTime, mass_initial, delta_mass); 

	% Running the IMM
	record_beta = [];
	record_theta = [];
	for i = 1:timeLength
		[x_est_current, move_cov_est_current, parameters] = filterIMM{itrMissile}.runFilter(i+startTime-1,deltaTime,y_obs{itrMissile}(:,i+1),measure_cov,i);    
		record_beta(end+1) = parameters{1};
		record_theta(end+1) = parameters{2};
	end

	%% Backward filtering
    % Define final state
	x_final = filterIMM{itrMissile}.x_ests_combined{end};

    % Each filter represents UKF with candidate beta & candidate theta
	filters_backward = cell(1);
	filterIMM_backward = cell(numMissile, 1);
    
	for b_idx = 1:length(beta_possible)
		for t_idx = 1:length(theta_possible)
			filters_backward{3*(b_idx-1)+t_idx} = all_UKF_backward(F_inv,[beta_possible(b_idx), theta_possible(t_idx)],H,[DY],x_final,movement_cov_final,observationSize,alpha,beta,kappa,gamma,mass_final);
		end
	end

	% Generate IMM integrating all the UKFs (filters)
	filterIMM_backward{itrMissile} = all_IMM_backward(filters_backward, x_final, movement_cov_final, beta_possible, theta_possible, startTime, endTime, mass_final, delta_mass); 

	% Running the IMM
	record_beta_backward = [];
	record_theta_backward = [];
	for i = timeLength-1:-1:1
		[x_est_current, move_cov_est_current, parameters] = filterIMM_backward{itrMissile}.runFilter(i+startTime+1,-deltaTime,y_obs{itrMissile}(:,i+1),measure_cov, timeLength-1+1-i);    
		record_beta_backward(end+1) = parameters{1};
		record_theta_backward(end+1) = parameters{2};
	end

	%% Kalman Smoothing
	smoothed = zeros(DX, timeLength);
	P_smoothed = cell(1, timeLength);
	smoothed(:, end) = filterIMM{itrMissile}.x_ests_combined{end};
	P_smoothed{end} = filterIMM{itrMissile}.move_cov_ests_combined{end};

	tempT = (timeLength-1) + 1;
	for i = timeLength-1:-1:1
		temp_z = inv(filterIMM_backward{itrMissile}.move_cov_ests_combined_prior{tempT-i}) * filterIMM_backward{itrMissile}.x_ests_combined_prior{tempT-i};
		temp_S = inv(filterIMM_backward{itrMissile}.move_cov_ests_combined_prior{tempT-i});
		Smoother_gain = filterIMM{itrMissile}.move_cov_ests_combined{i} * temp_S * inv(eye(DX) + filterIMM{itrMissile}.move_cov_ests_combined{i}*temp_S);
		P_smoothed{i} = filterIMM{itrMissile}.move_cov_ests_combined{i} - filterIMM{itrMissile}.move_cov_ests_combined{i} * temp_S * inv(eye(DX) + filterIMM{itrMissile}.move_cov_ests_combined{i}*temp_S) * filterIMM{itrMissile}.move_cov_ests_combined{i};
		smoothed(:,i) = filterIMM{itrMissile}.x_ests_combined{i} + (P_smoothed{i}*temp_z - Smoother_gain*filterIMM{itrMissile}.x_ests_combined{i});
	end

	%% Smoothing Likelihood
	obs_enu = zeros(DY, timeLength);

    for i = 1:timeLength
        obs_enu(:,i) = feval(H_inv, y_obs{1}(:,i+1));
    end
    
	smoothing_likelihood = 0;
	for i = 1:timeLength
		temp_likelihood = cal_smoothing_likelihood(obs_enu(:,i), smoothed(:,i), P_smoothed{i}, DY);
		smoothing_likelihood = smoothing_likelihood + temp_likelihood;
    end
    smoothing_likelihood = -smoothing_likelihood;

end