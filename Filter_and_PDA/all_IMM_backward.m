classdef all_IMM_backward < handle
    
    properties
        % properties of filter
        D                       % State dimension
        M                       % Measure dimension
        x_est                   % Estimated State
        move_cov_est            % Process Covariance
    
        % properties of IMM
        filters                 % Filters (UKFs)
        numFilters              % number of filters in IMM
        mu                      % cell of the mode probability of each time
        T
        x_ests                  % (cell of) Estimated States
        move_cov_ests           % (cell of) Process Covariacnes
        randas                  % (used for each filter)
        betas                   % (used for each filter)
        thetas                  % (used for each filter)
        residuals
        
        startTime               % Time when algorithm starts
        endTime                 % Time when algorithm ends

        massFinal               % Input final_mass value
		deltaMass               % Input delta_mass value
        
        % by YOON
        x_ests_prior
        move_cov_ests_prior
        x_ests_combined_prior
        move_cov_ests_combined_prior
        x_ests_combined         % (cell of) combined Estimated States
        move_cov_ests_combined  % (cell of) combined Estimated States
        % by YOON
    end
    
%% functions of IMM
    methods
        %% func1. Initialization
        function this = all_IMM_backward(filters, initialState, move_cov_est_initial, betas, thetas, startTime, endTime, mass_final, delta_mass)
            this.filters = filters;
            this.betas = betas;
            this.thetas = thetas;
            this.numFilters = length(filters);
            this.D = size(initialState,1);
            this.x_est = initialState;
            this.move_cov_est = move_cov_est_initial;
            this.startTime = startTime;
            this.endTime = endTime;
            
            if nargin == 9
                this.massFinal = mass_final;
				this.deltaMass = delta_mass;
            else
                fprintf('Oooooooops!');
                this.mass_cheating = 'none';
            end

            this.mu = cell(1);
            this.mu{1} = ones(this.numFilters,1) .* 1/this.numFilters;
            
            this.x_ests = cell(this.numFilters,1);
            this.move_cov_ests = cell(this.numFilters,1);
            for i = 1:this.numFilters
                this.x_ests{i,1} = initialState;
                this.move_cov_ests{i,1} = move_cov_est_initial;
            end
            
            % by YOON
            this.x_ests_prior = cell(this.numFilters, 1);
            this.move_cov_ests_prior = cell(this.numFilters, 1);
            for i = 1:this.numFilters
                this.x_ests_prior{i,1} = initialState;
                this.move_cov_ests{i,1} = move_cov_est_initial;
            end
            this.x_ests_combined_prior{1} = initialState;
            this.move_cov_ests_combined_prior{1} = move_cov_est_initial;
            this.x_ests_combined{1} = initialState;
            this.move_cov_ests_combined{1} = move_cov_est_initial;
            % by YOON
        end
        
        %% func2. Run Filter
        function [x_est_current,move_cov_est_current,params] = runFilter(this,time,deltaTime,y_obs,measure_cov,itr,...
                x_est_imm,move_cov_est_imm,probAssociation,probAssociationZero)
            
            % Our code : nargin = 5
            if nargin == 8
                this.x_est = x_est_imm;
                this.move_cov_est = move_cov_est_imm;
                fprintf('oops!!!!!!!!!!!!!!!!')
            end
            
            this.M = size(y_obs,1);
%             fprintf('time:%d -> %d\n',time,time+deltaTime*10);
            
          %% parameter setting for base filters
            % Beta setting for base filters

            mass = this.massFinal + (itr-1)*this.deltaMass;
            betas = this.betas;
            thetas = this.thetas;
            
            % Fill the params of each UKFs (filters)
            for b_idx = 1:length(this.betas)
                for t_idx = 1:length(this.thetas)
                    this.filters{3*(b_idx-1)+t_idx}.paramF(1) = betas(b_idx);
                    this.filters{3*(b_idx-1)+t_idx}.paramF(2) = thetas(t_idx);
                    this.filters{3*(b_idx-1)+t_idx}.mass = mass;
                end
            end
                      
          %% STEP 1 : Calculate the mixing probabilities
            % 1) T_k(i,j) = transition probability from mode_i to mode_j
            T_k = 50 * eye(this.numFilters) + ones(this.numFilters,this.numFilters);

            for i = 1:this.numFilters
                normalize = sum(T_k(i,:));
                T_k(i,:) = T_k(i,:) / normalize;
            end
            
            % 2) mu_k(i) = probability of mode_i
            k = length(this.mu);
            mu_k = this.mu{k};
            
            % 3) mu_small(i,j) = probability of being mode i now, then transitioning to mode_j
            mu_small = zeros(this.numFilters,this.numFilters);
            
            for i = 1:this.numFilters
                constant = 0;
                for j = 1:this.numFilters
                    % mu_small = TransitionProbability x mu_k
                    mu_small(i,j) = T_k(i,j)*mu_k(i);
                    constant = constant + mu_small(i,j);
                end
                for j = 1:this.numFilters
                    mu_small(i,j) = mu_small(i,j)/constant;
                end
            end
            
          %% STEP 2 : Calculate the mixed initial condition
            % mixed initial condition for the filter j at time k
            % 1) x_mixed : late states per a base filter adjusted by the current estimates (x_ests) and transitioning prob (mu_small)           
            x_mixed = cell(this.numFilters);
            for i = 1:this.numFilters
                x_mixed{i} = zeros(this.D,1);
                for j = 1:this.numFilters
                    x_mixed{i} = x_mixed{i} + this.x_ests{j,k} * mu_small(j,i);
                end
            end
            
            % 2) p_mixed : state covariance adjustment same as x_mized in the above
            p_mixed = cell(this.numFilters);
            for i = 1:this.numFilters
                p_mixed{i} = zeros(this.D,this.D);
                for j = 1:this.numFilters
                    p_mixed{i} = p_mixed{i} + mu_small(j,i)*(this.move_cov_ests{j,k}+(this.x_ests{j,k}-x_mixed{i})*(this.x_ests{j,k}-x_mixed{i})');
                end
            end
            
           %% STEP 3 : Perform mode_matched filtering
            for i = 1:this.numFilters
                [x_est_single_prior, move_cov_est_single_prior, x_est_single,move_cov_est_single,randa,residual] = this.filters{i}.runFilter(...
                    time,deltaTime,y_obs,measure_cov,x_mixed{i},p_mixed{i});
                this.x_ests{i,k+1} = x_est_single;
                this.move_cov_ests{i,k+1} = move_cov_est_single;
                this.randas(i,k+1) = randa;
                this.residuals(i,k+1) = residual;
                
                %by YOON
                this.x_ests_prior{i,k+1} = x_est_single_prior;
                this.move_cov_ests_prior{i, k+1} = move_cov_est_single_prior;
                % by YOON
            end
            
           %% STEP 4 : Update mode probability
            % 1) mu_k_1 : next time mode probability weighted by Randa
            mu_k_1 = zeros(this.numFilters,1);
            for i = 1:this.numFilters
                mu_k_1(i) = 0;
                for j = 1:this.numFilters
                    mu_k_1(i) = mu_k_1(i) + T_k(j,i)*mu_k(j);
                end
                mu_k_1(i) = mu_k_1(i) * this.randas(i,k+1);
            end
            mu_k_1 = mu_k_1 / sum(mu_k_1(:));
            
            % 2) expand mu by adding mode probability at tiem k+1
            this.mu{k+1} = mu_k_1;
            
           %% STEP 5 : Combine mode-conditioned estimates and covariances
            % 1) Latent state estimation at the IMM level.
            %    Weighted average of latent state sets from base filters.
            %    Weight is the mode probability.            
            this.x_est = zeros(this.D,1);
            for i = 1:this.numFilters
                this.x_est = this.x_est + mu_k_1(i)*this.x_ests{i,k+1};
            end

            % 2) Latent state covariance estimation at the IMM level.
            %    Weighted average of latent state sets from base filters.
            %    Weight is the mode probability.
            this.move_cov_est = zeros(this.D,this.D);
            for i = 1:this.numFilters
                this.move_cov_est = this.move_cov_est+mu_k_1(i)*(this.move_cov_ests{i,k+1}+(this.x_ests{i,k}-this.x_est)*(this.x_ests{i,k}-this.x_est)');
            end            

            % by YOON
            x_prior_temp = zeros(this.D,1);
            for i = 1:this.numFilters
                x_prior_temp = x_prior_temp + this.mu{k}(i)*this.x_ests_prior{i,k+1};
            end
            move_cov_prior_temp = zeros(this.D,this.D);
            for i = 1:this.numFilters
                move_cov_prior_temp = move_cov_prior_temp+this.mu{k}(i)*(this.move_cov_ests_prior{i,k+1}+(this.x_ests_prior{i,k}-x_prior_temp)*(this.x_ests_prior{i,k}-x_prior_temp)');
            end
            this.x_ests_combined_prior{k+1} = x_prior_temp;
            this.move_cov_ests_combined_prior{k+1} = move_cov_prior_temp;
            this.x_ests_combined{k+1} = this.x_est;
            this.move_cov_ests_combined{k+1} = this.move_cov_est;
            % by YOON
            
          %% OUTPUT : Return Purpose
            x_est_current = this.x_est;
            move_cov_est_current = this.move_cov_est;
            
            params = {};
            total_betas = reshape(cat(1, betas, betas, betas), 1, this.numFilters) * mu_k_1;
            total_thetas = cat(2, thetas, thetas, thetas) * mu_k_1;
            params{1} = total_betas;
            params{2} = total_thetas;
            
        end
    end
end
