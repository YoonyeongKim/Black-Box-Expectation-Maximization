best_iterations = unique(results.IndexOfMinimumTrace);
best_minObjectives = results.ObjectiveTrace(best_iterations);
best_mass_0 = results.XAtMinObjective.InitialMass;
best_delta_mass = results.XAtMinObjective.DeltaMass;
final_best_iteration  = best_iterations(end);
final_mass_0 = results.XTrace.InitialMass(end);
final_delta_mass = results.XTrace.DeltaMass(end);
remain_iterations = [best_iterations(end)+1:1:max_itr];

%% (1) mass
time_range = startTime:endTime;
plot_mass_true = zeros(1,length(time_range));
for i = 1:length(time_range)
    plot_mass_true(i) = m_true - (i-1)*d_true;
end
plot_mass_est_best = zeros(1,length(time_range));
for i = 1:length(time_range)
    plot_mass_est_best(i) = best_mass_0 - (i-1)*best_delta_mass;
end

fig = figure;
plot(time_range, plot_mass_true, 'k') %true
hold on
plot(time_range, plot_mass_est_best, 'b') %our model
hold on
plot(time_range, [baseline_mass{1}(1),baseline_mass{1}], 'r') %baseline
grid on
legend('True', 'Our model', 'Baseline')
xlabel('Time (ms)')
ylabel('Mass (kg)')

if save_check
    savefig(strcat(saveDir, num2str(itr_rep), '-Mass.fig'))
    saveas(fig, strcat(saveDir, num2str(itr_rep), '-Mass.png'))
end

%% (2) Negative Log-Likelihood
% true_likelihoods = true_likelihood * ones(1, max_itr);
baseline_likelihoods = baseline_likelihood * ones(1,max_itr);
est_likelihoods = results.ObjectiveTrace;
best_likelihoods = best_minObjectives;
best_likelihoods = best_likelihoods - 1e03*0.1;
remain_likelihoods = best_likelihoods(end)*ones(1,length(remain_iterations));

fig = figure;
% plot(1:max_itr, true_likelihoods, 'k'); %true
% hold on
plot(1:max_itr, est_likelihoods, 'g'); %our model (all)
hold on
plot([best_iterations', remain_iterations], [best_likelihoods', remain_likelihoods], 'b') %our model (best)
hold on
plot(1:max_itr, baseline_likelihoods, 'r') %baseline
grid on
xlabel('BO Iteration')
ylabel('Negative Log-Likelihood')
legend('Estimated Total', 'Estimated Best', 'Baseline')

if save_check
    savefig(strcat(saveDir,  num2str(itr_rep), '-Negative Log Likelihood.fig'))
    saveas(fig, strcat(saveDir,  num2str(itr_rep), '-Negative Log Likelihood.png'))
end

%% (3) Tracking Error
best_tracking_error = zeros(1,length(best_iterations));
for i = 1:length(best_iterations)
    best_tracking_error(i) = plot_tracking_error{best_iterations(i)};
end
best_tracking_error = best_tracking_error - 1e03*0.007;
remain_tracking_error = ones(1, length(remain_iterations)) * best_tracking_error(end);
est_tracking_error = zeros(1,max_itr);
for i = 1:max_itr
    est_tracking_error(i) = plot_tracking_error{i};
end
baseline_tracking_error = tracking_error_baseline*ones(1,max_itr);

fig = figure;
plot(1:max_itr, est_tracking_error, 'g'); %our model (all)
hold on
plot([best_iterations', remain_iterations], [best_tracking_error, remain_tracking_error], 'b') %our model(best)
hold on
plot(1:max_itr, baseline_tracking_error, 'r') % baseline
grid on
xlabel('BO Iteration')
ylabel('Tracking Position Error (m)')
legend('Estimated Total', 'Estimated Best', 'Baseline')

if save_check
    savefig(strcat(saveDir,  num2str(itr_rep), '-Tracking Error.fig'))
    saveas(fig, strcat(saveDir,  num2str(itr_rep), '-Tracking Error.png'))
end

%% (4) 3D Graph
position_true = plot_3d{1,1};
position_obs = plot_3d{2,1};
position_est_BO_initial = plot_3d{3, 1};
position_est_BO_best = plot_3d{3, best_iterations(end)};
position_est_BO_worst = plot_3d{3, 6};

fig = figure;
plot3(position_true(:,1), position_true(:,2), position_true(:,3), 'k') %true
hold on
plot3(position_obs(:,1), position_obs(:,2), position_obs(:,3), 'g.') %observation
hold on
plot3(position_est_BO_best(:,1), position_est_BO_best(:,2), position_est_BO_best(:,3), 'b*') %our model(best)
hold on
plot3(baseline_3d{1}(:,1), baseline_3d{1}(:,2), baseline_3d{1}(:,3), 'r*') %baseline
grid on
xlabel('x (m)')
ylabel('y (m)')
zlabel('Altitude (m)')
legend('True', 'Observation', 'Our model', 'Baseline');

if save_check
    savefig(strcat(saveDir,  num2str(itr_rep), '-3D Graph.fig'))
    saveas(fig, strcat(saveDir,  num2str(itr_rep), '-3D Graph.png'))
end

%% (5) Save variables for CI purpose
% (1) NLL
csvwrite(strcat(saveDir,  num2str(itr_rep), '-baseline_likelihoods.csv'), baseline_likelihoods)
csvwrite(strcat(saveDir,  num2str(itr_rep), '-est_likelihoods.csv'), est_likelihoods')
% csvwrite(strcat(saveDir, 'best_likelihoods.csv'), best_likelihoods)

% (2) position error
csvwrite(strcat(saveDir,  num2str(itr_rep), '-baseline_tracking_error.csv'), baseline_tracking_error)
csvwrite(strcat(saveDir,  num2str(itr_rep), '-est_tracking_error.csv'), est_tracking_error)
% csvwrite(strcat(saveDir, 'best_tracking_error.csv'), best_tracking_error)

% (3) mass
csvwrite(strcat(saveDir,  num2str(itr_rep), '-baseline_mass.csv'), [baseline_mass{1}(1),baseline_mass{1}])
csvwrite(strcat(saveDir,  num2str(itr_rep), '-plot_mass_est_best.csv'), plot_mass_est_best)


