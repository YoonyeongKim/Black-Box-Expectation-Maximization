classdef all_UKF_forward < handle
    
    properties
        % properties of filter
        D                       % State dimension
        M                       % Measure dimension
        F                       % Dynamics function
        paramF                  % Dynamics function parameter
        H                       % Measurement function
        paramH                  % Measurement function parameter
        x_est                   % State estimation
        move_cov_est            % process covariance
    
        % properties of UKF
        observationSize         % number of observations
        alpha                   % (related to spreas of the sample points around the mean)
        beta                    % (incorporates knowledge of prior distributions)
        kappa                   % (usually zero)
        randa                   % scaling parameter (decided by alpha, kappa, and observationSize)
        gamma                   % 
        etta_m                  % weight for mean estimation
        etta_c                  % weight for covariance estimation
        mass

    end
    
%% functions of UKF
    methods
       %% func1. Initialization
        function this = all_UKF_forward(F,paramF,H,paramH,initialState,move_cov_est_initial,observationSize,alpha,beta,kappa,gamma,mass)
            % initial setting for parameters
            this.D = size(initialState,1);
            this.x_est = initialState;
            this.F = F;
            this.paramF = paramF;
            this.H = H;
            this.paramH = paramH;
            this.move_cov_est = move_cov_est_initial;
            
            this.observationSize = observationSize;
            
            this.alpha = alpha;
            this.beta = beta;
            this.kappa = kappa;
            this.gamma = gamma;
            this.mass = mass;
            
            this.randa = alpha^2*(this.D+kappa)-this.D;
            
            % calculate weight for every (2D+1) sigma points
            for i = 1:2*this.D+1
                if i == 1
                    this.etta_m(i) = this.randa/( this.D + this.randa );
                    this.etta_c(i) = this.randa/( this.D + this.randa ) + 1 - this.alpha^2 + this.beta;
                else
                    this.etta_m(i) = 1 / ( 2 * ( this.D + this.randa ) );
                    this.etta_c(i) = 1 / ( 2 * ( this.D + this.randa ) );
                end
            end

        end
        
       %% func2. Run Filter : main function (Prediction step + Updating step)
        function [x_propagated, move_cov_propagated, x_est_current,move_cov_est_current,randa,residual] = runFilter(this,time,deltaTime,y_obs,measure_cov,x_est_imm,move_cov_est_imm)
            if nargin == 7
                this.x_est = x_est_imm;
                this.move_cov_est = move_cov_est_imm;
            end
            [x_propagated,move_cov_propagated,y_propagated,measure_cov_est,cross_cov_est] = runFilterPrediction(this,time,deltaTime,measure_cov,x_est_imm,move_cov_est_imm);
            [x_est_current,move_cov_est_current,randa,residual] = runFilterUpdate(this,time,deltaTime,y_obs,measure_cov,x_propagated,move_cov_propagated,y_propagated,measure_cov_est,cross_cov_est);
        end
        
       %% func3. Prediction step
        function [x_propagated,move_cov_propagated,y_propagated,measure_cov_est,cross_cov_est] = runFilterPrediction(this,time,deltaTime,measure_cov,x_est_imm,move_cov_est_imm)
            if nargin == 6
                this.x_est = x_est_imm;
                this.move_cov_est = move_cov_est_imm;
            end

           %% STEP 1 : Sampling
            this.move_cov_est = validateCovMatrix(this.move_cov_est);

            sP = this.gamma*chol(this.move_cov_est, 'lower');
            
            x_samplePoints = [this.x_est, this.x_est*ones(1,this.D)+sqrt(this.D+this.randa)*sP, this.x_est*ones(1,this.D)-sqrt(this.D+this.randa)*sP];
                             % one point from mean
                             % D points from mean + chol term
                             % D points from mean - chol term
            x_samplePoints = real(x_samplePoints);

          %% STEP 2 : Prediction
            % propagated the sigma points through the prediction function F
            % 1) mean of the prediction
            x_samplePoints_propagated = zeros(this.D,2*this.D+1);
            for j = 1:(2*this.D+1)
                x_samplePoints_propagated(:,j) = feval(this.F, x_samplePoints(:,j), time*0.1, deltaTime, this.paramF(1), this.paramF(2), this.mass);
            end
            x_propagated = x_samplePoints_propagated*this.etta_m';

            % 2) covariance of the prediction
            move_cov_propagated = zeros(this.D,this.D);
            for j = 1:(2*this.D+1)
                move_cov_propagated = move_cov_propagated + this.etta_c(j) * (x_samplePoints_propagated(:,j)-x_propagated)*(x_samplePoints_propagated(:,j)-x_propagated)';
            end
            
           %% STEP 3 : Observation
            % propagated the sigma points through the observation function H
            % 1) mean of the observation
            y_samplePoints_propagated = zeros(this.observationSize,2*this.D+1);
            for j = 1:(2*this.D+1)
                y_samplePoints_propagated(:,j) = feval(this.H, x_samplePoints_propagated(:,j), -1, -1, this.paramH);
            end            
            y_propagated = y_samplePoints_propagated*this.etta_m';

            % 2) covariance of the observation
            measure_cov_est = measure_cov;
            for j = 1:(2*this.D+1)
                measure_cov_est = measure_cov_est + this.etta_c(j)*(y_samplePoints_propagated(:,j)-y_propagated(:))*(y_samplePoints_propagated(:,j)-y_propagated(:))';
            end
            
            % 3) cross covariance of the predicted observation and predicted state
            cross_cov_est = zeros(this.D,this.observationSize);
            for j = 1:(2*this.D+1)
                cross_cov_est = cross_cov_est + this.etta_c(j)*(x_samplePoints_propagated(:,j)-x_propagated(:))*(y_samplePoints_propagated(:,j)-y_propagated(:))';
            end            
        end
        
       %% func4. Updating step
        function [x_est_current,move_cov_est_current,randa,residual] = runFilterUpdate(this,time,deltaTime,y_obs,measure_cov,x_propagated,move_cov_propagated,y_propagated,measure_cov_est,cross_cov_est)
            this.M = size(y_obs,1);

            % compute Kalman Gain
            gain = cross_cov_est / measure_cov_est;
            
            % combine information above to find the estimation of mean and covariance
            this.x_est = x_propagated + gain*(y_obs-y_propagated);
            this.move_cov_est = move_cov_propagated - gain*measure_cov_est*gain';
            
            x_est_current = this.x_est;
            move_cov_est_current = this.move_cov_est;
            measure_cov_est(isnan(measure_cov_est))=0;
            measure_cov_est = measure_cov_est + 1*eye(this.M);
            randa = mvnpdf(y_obs-y_propagated,zeros(this.M,1),measure_cov_est);
            residual = sqrt( (y_obs-y_propagated)'*(y_obs-y_propagated) );
        end        
    end
    
end

