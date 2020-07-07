%% Readme
% This code is written to run Bayesian Optimization test to find the initial mass and the change rate of the mass.
% You can change the settings of this experiment as you want.

% 01. Directory for saving result figures (Line 23~27)
% 02. Data properties (Line 51~52, The variable name should be 'ballistic_target' to run properly.)
% 03. Index of vehicle (Line 54)
% 04. Time interval of the trajectory (Line 63~64)
% 05. Error level (Line 69~73)
% 06. Radar position (Line 76~78)
% 07. Parameters of UKF (Line 82~86)
% 08. Hypothesis value for mass, beta, and theta for BayesOpt Test (Line 99~111)
% 09. Covariance matrix of measurement in IMM (Line 138)
% 10. Covariance matrix of movement in IMM (Line 139~140)
% 11. Iteration for Bayesian Optimization test (Line 116)

%% Path Adjusting
clear
clc
cwd = pwd;
addpath(genpath(cwd));

save_check = true;
saveDir = 'Results\';
if save_check
    mkdir(saveDir);
end

%% Global Variables
% missile
global idxMissile numMissile curMissile
% time
global startTime endTime timeLength deltaTime
% error level
global startErrorLevel
% experiment setting
global DX DY beta_possible theta_possible mass_possible y_obs y_trues data_name ballistic_data
% IMM matrix
global measure_cov movement_cov_initial movement_cov_final
% radar position
global rx ry rz
% UKF and IMM parameter
global alpha beta kappa gamma observationSize itrMissile
% plot purpose
global BO_itr plot_tracking_error plot_IMM plot_enu plot_3d
% baseline purpose
global baseline_likelihood baseline_mass tracking_error_baseline baseline_3d

%% Settings for Experiment
% Load data
data_name = 'ballistic_target_YY_180518.mat';
load(data_name);

idxMissile = [1];
numMissile = length(idxMissile);
curMissile = 1;

ballistic_data = ballistic_target(idxMissile(curMissile));
DX = 9;
DY = 3;

% Define time interval
startTime = 410;
timeLength = 100;
endTime = startTime + timeLength;
deltaTime = 0.1;

% Define noise factors
noise_level = 0.05;
startErrorLevel = 0.003;
angleNoiseLevel = 0.0001*pi*2;
positionNoiseLevel = 1;
noise_constant = 1000;

% Define radar position
rx = -134440;
ry = -19870;
rz = -1410;

%% Settings for UKF and IMM
% UKF
alpha = 0.15;            % distance from the centroid to the sigma points and the weights of the sigma points
beta = 2;                % the distribution function ( 2 == Gaussian distribution)
kappa = 0;               % usually set to zero
gamma = 1;               % the spread of the the sigma points
observationSize = DY;    % the dimension of the obesrvations

% IMM
itrMissile = curMissile;

%% Settings for Data

% Make true trajectory data
trajectories_true{1} = ballistic_data.x;
trajectories_true_cut{1} = trajectories_true{1}(:, startTime:endTime);
T = size(trajectories_true_cut{1},2);

% Make interval of mass, beta, and theta for BayesOpt Test
m_true = ballistic_data.total_mass(startTime);
m_lb = m_true - 0.05e3;
m_ub = m_true + 0.05e3;

d_true = (ballistic_data.total_mass(startTime)-ballistic_data.total_mass(endTime))/timeLength;
d_lb = d_true - 1;
d_ub = d_true + 1;

beta_true = mean(ballistic_data.ballistic_coeff(startTime:endTime));
theta_true = mean(ballistic_data.magThrBMArr(startTime:endTime));
beta_possible = [beta_true-0.5*1e04, beta_true, beta_true+0.5*1e04];
theta_possible = [theta_true-0.003*1e05, theta_true, theta_true+0.003*1e05];
mass_possible = [m_true-0.05e3, m_true, m_true+0.05e3];

%% Repeated Experiment

num_rep = 5;
max_itr = 10;
save_baseline_likelihoods = zeros(num_rep, max_itr);
save_est_likelihoods = zeros(num_rep, max_itr);
save_baseline_tracking_error = zeros(num_rep, max_itr);
save_est_tracking_error = zeros(num_rep, max_itr);
save_baseline_mass = zeros(num_rep, timeLength+1);
save_est_mass = zeros(num_rep, timeLength+1);

for itr_rep = 1:num_rep    
    randomseed = itr_rep;
    rand('seed', randomseed);
    randn('seed', randomseed);
    
    %% Global Variables
    baseline_mass = {};
    baseline_likelihood = 0;
    tracking_error_baseline = 0;
    baseline_3d = {};
    
    %% Run true_IMM for baseline
    % Settings for IMM
    noise_cov = [positionNoiseLevel 0 0; 0 angleNoiseLevel 0; 0 0 angleNoiseLevel];
    measure_cov = noise_cov + noise_constant * ones(DY,DY);
    movement_cov_initial = 10000*eye(DX)+1000*ones(DX,DX);
    movement_cov_final = movement_cov_initial;

    % Prepare true and observation data by the radar
    y_trues = cell(numMissile,1);
    y_obs = cell(numMissile,1);
    for k = 1:numMissile
        y_trues{k} = zeros(DY,T);
        y_obs{k} = zeros(DY,T);
        for i = 1:T
            obs = dyn_ballisticM_HV61_Measurement_YY(trajectories_true_cut{k}(:,i),-1,-1,DY);
            y_trues{k}(:,i) = obs;
            y_obs{k}(:,i) = addRadarNoise(y_trues{k}(:,i),noise_level);
        end
    end

    %% Run true_IMM
    true_likelihood = runKalmanSmoothing_for_true_mass(m_true, d_true);

    %% Run baseline
    temp = runKalmanSmoothing_baseline;

    %% Run BO Test    

    % For plot purpose
    BO_itr = 0;
    plot_tracking_error = {};
    plot_IMM = {};
    plot_enu = {};
    plot_3d = {};

    % Run Bayesian Optimization Test
    var_mass_0 = optimizableVariable('InitialMass', [m_lb, m_ub], 'Type', 'real');
    var_delta_mass = optimizableVariable('DeltaMass', [d_lb, d_ub], 'Type', 'real');

    vars = [var_mass_0,var_delta_mass];
    fn = @(vars)runKalmanSmoothing(vars.InitialMass, vars.DeltaMass);
    results = bayesopt(fn, vars, 'Verbose',...
        1,'AcquisitionFunctionName','expected-improvement', 'MaxObjectiveEvaluations', max_itr)
    
    %% Save Result Figures
    draw_chart;    
    close all
    
end