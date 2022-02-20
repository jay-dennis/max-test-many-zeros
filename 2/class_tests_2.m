classdef class_tests_2 < class_dgp_2
    properties (SetAccess=public, GetAccess=public)
        % Max test and max t-test variables
        k_theta_n
        theta_star
        y_hat
        % V_all, V_theta
        max_test_stat, max_t_test_stat
        % T_max_test, T_max_t_test
        p_max_test, p_max_t_test
        p_max_test_bs, p_max_t_test_bs
        dr_max
        dr_max_t
        M_max_sim
        % Asymptotic Wald Test
        T_wald
        p_wald
        % Parametric Bootstrapped Wald Tests
        M
        % T_star_wald_bs
        p_wald_bs
        % other
        alpha_levels
        dr
        p_vals
    end
    properties (SetAccess=public, GetAccess=public, Hidden = true)
    end
    methods
        function out = about(obj)
            display('Generates 6.2: Tests of Linear Regression Zero Restrictions');
            display('Version 6');
            display('------------');
        end
        function obj = class_tests_2(dgp_type, hypothesis_type, n, seed_input, k_theta, k_delta)
            obj@class_dgp_2(dgp_type, hypothesis_type, n, seed_input, k_theta, k_delta);
        end
        function obj = set_params(obj, k_theta_n)
            obj.alpha_levels = [0.01 0.05 0.10];
            obj.k_theta_n = k_theta_n;
            obj.M_max_sim=25000; % 25000;
            obj.M = 1000;
        end
        % All tests
        function obj = run_all_tests_fcn(obj, k_theta_n)
            obj = obj.set_params(k_theta_n);
            obj = obj.max_test_stat_fcn;
            %
            [obj.T_wald, obj.p_wald] = obj.wald_fcn;
            [~, obj.p_wald_bs] = obj.bs_wald_fcn(obj.T_wald);
            %
            obj.p_vals = [obj.p_max_test; obj.p_max_t_test; obj.p_max_test_bs; obj.p_max_t_test_bs; obj.p_wald; obj.p_wald_bs];
            obj.dr = obj.set_decision_rules;
        end
        function obj = clean_up(obj)
            obj.Y0=[];
            obj.e=[];
            obj.x_delta=[];
            obj.x_theta=[];
            obj.y_hat=[];
            obj.theta_0 = [];
            obj.delta_0 = [];
            obj.theta_star=[];
        end
        function dr_vec = set_decision_rules(obj)
            p_temp = obj.p_vals;
            len   = length(p_temp);
            len_a = length(obj.alpha_levels);
            dr_vec = zeros(len,len_a);
            for a = 1:len_a
                alpha_level = obj.alpha_levels(a);
                for l = 1:len
                    p = p_temp(l);
                    if p >= alpha_level
                        dr_temp = 0; % Fail to reject H_0
                    elseif p < alpha_level
                        dr_temp = 1; % Reject H_0
                    end
                    dr_vec(l,a) = dr_temp;
                end
            end
        end
        
        % Max Test Methods
        function obj = max_test_stat_fcn(obj)
            obj.theta_star = zeros(obj.k_theta_n,(obj.k_delta + 1));
            obj.y_hat = zeros(obj.n,obj.k_theta_n);
            for i = 1:obj.k_theta_n
                temp_x = [obj.x_delta; obj.x_theta(i, :)];
                temp_theta = temp_x' \ obj.Y0'; % Just OLS here;
                obj.theta_star(i,:) = temp_theta;
                obj.y_hat(:,i) = temp_theta' * temp_x;
            end
            [V_all, V_theta] = obj.V_hat_fcn;
            
            W = ones(obj.k_theta_n,1);
            obj.max_test_stat = sqrt(obj.n) * max(abs(W .* obj.theta_star(:,end)));
            Wt = sqrt(diag(V_theta)).^(-1);
            obj.max_t_test_stat = sqrt(obj.n) * max(abs(Wt .* obj.theta_star(:,end)));
            
            W_in = [W Wt]; stat_in = [obj.max_test_stat obj.max_t_test_stat];
            [~, p_out] = obj.max_test_pval_fcn(V_all, obj.M_max_sim, W_in, stat_in);
            obj.p_max_test = p_out(1); obj.p_max_t_test = p_out(2);
            
            [~, p_out] = obj.max_test_bs_pval_fcn(obj.M_max_sim, W_in, stat_in);
            obj.p_max_test_bs = p_out(1); obj.p_max_t_test_bs = p_out(2);
        end
        function [V, V_theta] = V_hat_fcn(obj)
            V = zeros(obj.k_theta_n*(obj.k_delta + 1), obj.k_theta_n*(obj.k_delta + 1));
            V_theta = zeros(obj.k_theta_n,obj.k_theta_n);
            for i = 1:obj.k_theta_n
                x_i = [obj.x_delta; obj.x_theta(i, :)];
                delta_hat_i = obj.theta_star(i,1:end-1)';
                J_i = (x_i * x_i') / obj.n;
                temp_e_hat_i = obj.Y0 - delta_hat_i' * obj.x_delta;
                for j = 1:obj.k_theta_n
                    x_j = [obj.x_delta; obj.x_theta(j, :)];
                    delta_hat_j = obj.theta_star(j,1:end-1)';
                    J_j = (x_j * x_j') / obj.n;
                    temp_e_hat_j = obj.Y0 - delta_hat_j' * obj.x_delta;
                    temp_i = repmat(temp_e_hat_i, obj.k_delta+1, 1) .* x_i;
                    temp_j = repmat(temp_e_hat_j, obj.k_delta+1, 1) .* x_j;
                    S = temp_i * temp_j' / obj.n; 
                    v = ( J_i \ eye(size(J_i,1)) ) * S * ( J_j \ eye(size(J_j,1)) );
                    V(((i-1)*(obj.k_delta + 1)+1):(i*(obj.k_delta + 1)),((j-1)*(obj.k_delta + 1)+1):(j*(obj.k_delta + 1))) = v;
                    V_theta(i,j) = v(end,end);
                end
            end
        end
        function [T_out, p_out] = max_test_pval_fcn(obj, V, M, W_in, stat_in)
            %temp_v = chol(V);
            temp_v = sqrtm(V);
            Z = randn(obj.k_theta_n*(obj.k_delta + 1),M);
            Z_sim = temp_v * Z;
            Z_theta = zeros(obj.k_theta_n,M);
            for i = 1:obj.k_theta_n
                Z_theta(i,:) = Z_sim(i*(obj.k_delta + 1),:);
            end
            T_out = zeros(length(stat_in), M); p_out = zeros(length(stat_in),1);
            for j = 1:length(stat_in)
                W = W_in(:, j); statistic = stat_in(j);
                W = repmat(W, 1, M);
                T = max(abs(W .* Z_theta),[],1);
                p = obj.pval_fcn(T, M, statistic);
                T_out(j,:) = T; p_out(j) = p;
            end
        end
        function [T_out, p_out] = max_test_bs_pval_fcn(obj, M, W_in, stat_in)
            M = 1000;
            % Bootstrapped Max Test
            % estimate imposing the null.
            temp_X = obj.x_delta';
            % d_hat_null = inv(temp_X' * temp_X) * temp_X' * obj.Y0';
            d_hat_null = temp_X \ obj.Y0'; % faster
            % form null imposed residuals.
            y_hat_null = d_hat_null' * obj.x_delta;
            e_hat_null = obj.Y0 - y_hat_null;
            % mutliply resids by normal rv
            Z = randn(M, obj.n);
            e_hat_wbs = repmat(e_hat_null, M, 1) .* Z;
            % construct y*
            y_wbs = repmat(y_hat_null, M, 1) + e_hat_wbs;
            % estimate pars models
            theta_hat_wbs = zeros(M, obj.k_theta_n);
            for m = 1:M
                temp_Y = y_wbs(m,:)';
                for i = 1:obj.k_theta_n
                    temp_X = [obj.x_delta' obj.x_theta(i,:)'];
                    % b_hat_wbs = inv(temp_X' * temp_X) * temp_X' * temp_Y;
                    b_hat_wbs = temp_X \ temp_Y;
                    theta_hat_wbs(m,i) = b_hat_wbs(end);
                end
            end
            % form max stats, and sort
            % feed to pval func
            T_out = zeros(M, length(stat_in)); p_out = zeros(length(stat_in),1);
            for j = 1:length(stat_in)
                W = W_in(:, j); statistic = stat_in(j);
                W = repmat(W', M, 1);
                T = sqrt(obj.n) * max(abs(W .* theta_hat_wbs),[],2);
                T = sort(T);
                p = obj.pval_fcn(T, M, statistic);
                T_out(:,j) = T; p_out(j) = p;
            end
        end
       
        % Wald
        function [wald, p_val] = wald_fcn(obj)
            y0 = obj.Y0';
            X = [obj.x_delta' obj.x_theta(1:obj.k_theta_n,:)'];
            temp_beta = X \ y0; % Just OLS here;
            e_hat = y0 -  X * temp_beta;
            % Use White's Robust variance estimator 
            temp = X .* repmat(e_hat, 1, obj.k_delta + obj.k_theta_n);
            % var_B_White = inv(X' * X) * (Z' * Z) * inv(X' * X);
            outer = ((X' * X) \ eye(obj.k_delta + obj.k_theta_n));
            inner = temp' * temp;
            v = outer * inner * outer;
            v_theta = v((obj.k_delta+1):end,(obj.k_delta+1):end);
            v_inv = (v_theta \ eye(obj.k_theta_n));
            %
            theta_vec = temp_beta((obj.k_delta+1):end);
            wald = (theta_vec)' * v_inv * (theta_vec);
            p_val = 1 - chi2cdf(wald, obj.k_theta_n);
        end
        % Bootstrapped Wald
        function [T_star, p_val] = bs_wald_fcn(obj, wald)
            % impose the null first
            X_null = obj.x_delta'; X_all = [obj.x_delta' obj.x_theta(1:obj.k_theta_n,:)'];
            y0 = obj.Y0';
            beta_hat_null = X_null \ y0; % Just OLS here;
            y_hat_null = X_null * beta_hat_null;
            e_null = y0 - y_hat_null;
            % generate null imposed y's
            Z = randn(obj.n, obj.M);
            e_bs = repmat(e_null, 1, obj.M) .* Z;
            y_null_bs = y_hat_null + e_bs;
            %
            T_star = zeros(obj.M,1);
            for m = 1:obj.M
                y_temp = y_null_bs(:,m);
                b_hat = X_all \ y_temp;
                theta_vec = b_hat(obj.k_delta+1:end);
                e_hat = y_temp - X_all * b_hat;
                %
                temp = X_all .* repmat(e_hat, 1, obj.k_delta + obj.k_theta_n);
                outer = ((X_all' * X_all) \ eye(obj.k_delta + obj.k_theta_n));
                inner = temp' * temp;
                v = outer * inner * outer;
                %
                v_theta = v((obj.k_delta+1):end,(obj.k_delta+1):end);
                v_inv = (v_theta \ eye(obj.k_theta_n));
                T_star(m) = (theta_vec)' * v_inv * (theta_vec);
            end
            T_star = sort(T_star);
            p_val = obj.pval_fcn(T_star, obj.M, wald);
        end
        
        % General Functions
        function p = pval_fcn(obj, T, M, statistic)
            p = (1/M) * sum((T>statistic));
        end
    end
end
