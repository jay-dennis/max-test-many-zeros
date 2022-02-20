classdef class_tests_1 < class_dgp_1
    properties (SetAccess=public, GetAccess=public)
        % Max test and max t-test variables
        k_theta_n
        k_theta_n_vec
        theta_star
        y_hat
        % V_theta_hat, V_theta_HAC_hat
        max_test_stat_vec, max_t_test_stat_vec, max_t_test_stat_HAC_vec
        p_max_test_vec, p_max_t_test_vec, p_max_t_test_HAC_vec
        % T_max_test, T_max_t_test, T_max_t_test_HAC
        p_max_test_bs_vec, p_max_t_test_bs_vec, p_max_t_test_HAC_bs_vec
        % T_max_test_bs, T_max_t_test_bs, T_max_t_test_HAC_bs
        dr_max
        dr_max_t
        dr_max_t_HAC
        M_max_sim
        % Ljung-Box
        LBQ_test_stat_vec
        % LBQ_DWB
        p_LBQ_vec
        dr_LBQ
        % AP sup LM
        lambda_grid
        sup_LM_test_stat
        % sup_LM_DWB
        p_sup_LM
        dr_sup_LM
        % CvM test
        pi_grid
        CvM_test_stat
        % CvM_DWB
        p_CvM
        dr_CvM
        % Max Correlation Test
        T_max_corr_vec
        % T_max_corr_DWB
        p_max_corr_vec
        dr_max_corr
        % DWB
        bn
        M
        % other
        alpha_levels
        dr, dr_vec
        p_vals, p_vals_vec
    end
    properties (SetAccess=public, GetAccess=public, Hidden = true)
    end
    methods
        function out = about(obj)
            display('Generates 6.1: White Noise Test');
            display('Version 6');
            display('------------');
        end
        function obj = class_tests_1(dgp_type, n, seed_input)
            obj@class_dgp_1(dgp_type, n, seed_input);
        end
        function obj = set_params(obj, k_theta_n_vec)
            obj.alpha_levels = [0.01 0.05 0.10];
            obj.k_theta_n_vec = k_theta_n_vec;
            obj.k_theta_n = max(obj.k_theta_n_vec); % reduces redunancy
            obj.lambda_grid = (-.8:.005:.8); % obj.lambda_grid = (-.8:.005:.8);
            obj.pi_grid = (0:0.001:3.141); % obj.pi_grid = (0:0.001:3.141);
            obj.bn = floor(sqrt(obj.n));
            obj.M_max_sim=25000; %25000;
            obj.M = 1000;
        end
        % All tests
        function obj = run_all_tests_fcn(obj, k_theta_n_vec)
            obj = obj.set_params(k_theta_n_vec);
            obj = obj.max_test_stat_fcn;
            %
            [gamma_vec, rho_vec] = obj.fcn_rho_vec;
            obj.LBQ_test_stat_vec = obj.LBQ_fcn(rho_vec);
            obj.T_max_corr_vec = obj.max_corr_fcn(rho_vec);
            obj.sup_LM_test_stat = obj.sup_LM_fcn(rho_vec);
            obj.CvM_test_stat = obj.CvM_fcn(gamma_vec);
            %
            LBQ_DWB_vec = zeros(obj.M,length(obj.k_theta_n_vec));
            T_max_corr_DWB_vec = zeros(obj.M,length(obj.k_theta_n_vec));
            sup_LM_DWB = zeros(obj.M,1);
            CvM_DWB = zeros(obj.M,1);
            for m = 1:obj.M
                w = obj.DWB_w_fcn;
                [gamma_vec, rho_vec] = obj.fcn_rho_vec(w);
                LBQ_DWB_vec(m,:) = obj.LBQ_fcn(rho_vec);
                T_max_corr_DWB_vec(m,:) = obj.max_corr_fcn(rho_vec);
                sup_LM_DWB(m) = obj.sup_LM_fcn(rho_vec);
                CvM_DWB(m) = obj.CvM_fcn(gamma_vec);
            end
            
            for ind = 1:length(obj.k_theta_n_vec)
                obj.p_LBQ_vec{ind} = obj.pval_fcn(LBQ_DWB_vec(:,ind), obj.M, obj.LBQ_test_stat_vec(:,ind));
                obj.p_max_corr_vec{ind} = obj.pval_fcn(T_max_corr_DWB_vec(:,ind), obj.M, obj.T_max_corr_vec(:,ind));
            end
                obj.p_sup_LM = obj.pval_fcn(sup_LM_DWB, obj.M, obj.sup_LM_test_stat);
                obj.p_CvM = obj.pval_fcn(CvM_DWB, obj.M, obj.CvM_test_stat);
            %
            for ind = 1:length(obj.k_theta_n_vec)
                obj.p_vals_vec{ind} = [obj.p_max_test_vec{ind}; obj.p_max_t_test_vec{ind}; obj.p_max_t_test_HAC_vec{ind}; obj.p_max_test_bs_vec{ind}; obj.p_max_t_test_bs_vec{ind}; obj.p_max_t_test_HAC_bs_vec{ind}; obj.p_LBQ_vec{ind}; obj.p_sup_LM; obj.p_CvM; obj.p_max_corr_vec{ind}];
            end
        end
        function obj = clean_up(obj)
            obj.Y0=[];
            obj.e=[];
            obj.y_hat=[];
            obj.theta_star=[];
            obj.lambda_grid=[];
            obj.pi_grid=[];
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
        function obj = fcn_final_output(obj, k_theta_n)
            ind = find(obj.k_theta_n_vec == k_theta_n, 1);
            if isempty(ind) == 1
                fprintf('\n Error - enter a valid k_theta_n \n')
            else
                obj.k_theta_n = k_theta_n; 
                obj.p_vals = obj.p_vals_vec{ind};
                obj.dr = obj.set_decision_rules;
            end
        end
        
        % Max Test 
        function obj = max_test_stat_fcn(obj)
            obj.theta_star = zeros(obj.k_theta_n,2);
            obj.y_hat = zeros(obj.n,obj.k_theta_n);
            for i = 1:obj.k_theta_n
                temp_y = obj.Y0(1+i:end);
                temp_x = [ones(length(temp_y),1) obj.Y0(1:end-i)];
                temp_theta = temp_x \ temp_y; % Just OLS here;
                obj.theta_star(i,:) = temp_theta;
                obj.y_hat(1+i:end,i) = temp_x * temp_theta;
            end
            [V_theta_hat, V_theta_HAC_hat, V_all, V_all_HAC] = obj.V_hat_fcn;
            W = ones(obj.k_theta_n,1);
            Wt = sqrt(diag(V_theta_hat)).^(-1);
            Wt_HAC = sqrt(diag(V_theta_HAC_hat)).^(-1);
            
            for ind = 1:length(obj.k_theta_n_vec)
                k = obj.k_theta_n_vec(ind);
                theta_temp = obj.theta_star(1:k,2);
                W_temp = W(1:k); Wt_temp = Wt(1:k); Wt_HAC_temp = Wt_HAC(1:k);
                obj.max_test_stat_vec{ind} = sqrt(obj.n) * max(abs(W_temp .* theta_temp));
                obj.max_t_test_stat_vec{ind} = sqrt(obj.n) * max(abs(Wt_temp .* theta_temp));
                obj.max_t_test_stat_HAC_vec{ind} = sqrt(obj.n) * max(abs(Wt_HAC_temp .* theta_temp));
            end
            
            [~, obj.p_max_test_vec] = obj.max_test_pval_fcn(V_all, obj.M_max_sim, W, obj.max_test_stat_vec);
            [~, obj.p_max_t_test_vec] = obj.max_test_pval_fcn(V_all, obj.M_max_sim, Wt, obj.max_t_test_stat_vec);
            [~, obj.p_max_t_test_HAC_vec] = obj.max_test_pval_fcn(V_all_HAC, obj.M_max_sim, Wt_HAC, obj.max_t_test_stat_HAC_vec);

            %Using both right now to see which works better...
            W_in = [W Wt Wt_HAC]; 
            for ind = 1:length(obj.k_theta_n_vec)
                stat_in{ind} = [obj.max_test_stat_vec{ind} obj.max_t_test_stat_vec{ind} obj.max_t_test_stat_HAC_vec{ind}];
            end
            [~, p_out] = obj.max_test_pval_bootstrap_fcn(obj.M_max_sim, W_in, stat_in);
            for ind = 1:length(obj.k_theta_n_vec)
                obj.p_max_test_bs_vec{ind} = p_out{ind}(1); 
                obj.p_max_t_test_bs_vec{ind} = p_out{ind}(2); 
                obj.p_max_t_test_HAC_bs_vec{ind} = p_out{ind}(3);
            end
        end
        function [v_theta, v_theta_HAC, v_all, v_all_HAC] = V_hat_fcn(obj)
            bandwidth = floor(3 * obj.n^((3/10) - .00001));
            v_theta = zeros(obj.k_theta_n,obj.k_theta_n);
            v_theta_HAC = zeros(obj.k_theta_n,obj.k_theta_n);
            v_all = zeros(2*obj.k_theta_n,2*obj.k_theta_n);
            v_all_HAC = zeros(2*obj.k_theta_n,2*obj.k_theta_n);
            mean_Y0 = sum(obj.Y0) / obj.n;
            method = 3;
            if method == 1
                S = zeros(obj.k_theta_n,obj.k_theta_n);
                S_HAC = zeros(obj.k_theta_n,obj.k_theta_n);
                for i = 1:obj.k_theta_n
                    for j = 1:obj.k_theta_n
                        for t = (max(i,j)+1):obj.n
                            S(i,j) = S(i,j) + (1/obj.n) * ((obj.Y0(t) - mean_Y0).^2 * obj.Y0(t-i) * obj.Y0(t-j));
                        end
                        for s = (i+1):obj.n
                        for t = (j+1):obj.n
                            K = max(0,1 - abs((s-t)/bandwidth));
                            S_HAC(i,j) = S_HAC(i,j) + (1/obj.n) * K * ((obj.Y0(t) - mean_Y0) * (obj.Y0(s) - mean_Y0) * obj.Y0(s-i) * obj.Y0(t-j));
                        end
                        end
                        v_theta(i,j) = S(i,j) / ((1/obj.n) * sum(obj.Y0 .^2))^(2);
                        v_theta_HAC(i,j) = S_HAC(i,j) / ((1/obj.n) * sum(obj.Y0 .^2))^(2);
                    end
                end
            elseif method == 2
                S = zeros(obj.k_theta_n,obj.k_theta_n);
                S_HAC = zeros(obj.k_theta_n,obj.k_theta_n);
                Y0_demeaned = obj.Y0 - mean_Y0;
                temp_Y = obj.Y0;
                K_mat0 = zeros(obj.n, obj.n);
                temp = (1:obj.n)' - 1;
                for t = 1:obj.n-1
                    K_mat0(t:obj.n,t) = temp(1:end-t+1);
                end
                K_mat0 = K_mat0 + K_mat0';                
                K_mat0 = max(0,1 - (K_mat0./bandwidth));
                for i = 1:obj.k_theta_n
                    for j = 1:obj.k_theta_n
                        diff = abs(j-i);
                        if j >= i
                            y = Y0_demeaned(j+1:end); 
                            wi = temp_Y(diff+1:end-i); wj = temp_Y(1:end-j);
                            K_mat = K_mat0(j+1:end, j+1:end);
                        elseif j < i
                            y = Y0_demeaned(i+1:end); 
                            wi = temp_Y(1:end-i); wj = temp_Y(diff+1:end-j);
                            K_mat = K_mat0(i+1:end, i+1:end);
                        end
                        Si = y .* wi; Sj = y .* wj;
                        S(i,j) = Si' * Sj / obj.n;
                        S_HACij = Si * Sj' .* K_mat;
                        S_HAC(i,j) = sum(sum(S_HACij)) / obj.n;
                        denom = ( sum(obj.Y0 .^2) / obj.n )^(2);
                        v_theta(i,j) = S(i,j) / denom;
                        v_theta_HAC(i,j) = S_HAC(i,j) / denom;
                    end
                end
            elseif method == 3
                Y0_demeaned = obj.Y0 - mean_Y0;
                temp_Y = obj.Y0;
                K_mat0 = zeros(obj.n, obj.n);
                temp = (1:obj.n)' - 1;
                for t = 1:obj.n-1
                    K_mat0(t:obj.n,t) = temp(1:end-t+1);
                end
                K_mat0 = K_mat0 + K_mat0';                
                K_mat0 = max(0,1 - (K_mat0./bandwidth));
                for i = 1:obj.k_theta_n
                    for j = 1:obj.k_theta_n
                        diff = abs(j-i);
                        if j >= i
                            y = Y0_demeaned(j+1:end); 
                            wi = temp_Y(diff+1:end-i); wj = temp_Y(1:end-j);
                            K_mat = K_mat0(j+1:end, j+1:end);
                        elseif j < i
                            y = Y0_demeaned(i+1:end); 
                            wi = temp_Y(1:end-i); wj = temp_Y(diff+1:end-j);
                            K_mat = K_mat0(i+1:end, i+1:end);
                        end
                        one_vec = ones(length(y), 1);
                        xi = [one_vec wi]; xj = [one_vec wj];
                        Si = y .* xi; Sj = y .* xj;
                        S = Si' * Sj / obj.n;
                        S_HAC = Si' * K_mat * Sj / obj.n;
                        J = [1 (sum(y)/obj.n); (sum(y)/obj.n) (sum(y.^2)/obj.n)];
                        J_inv = J \ eye(2);
                        v = J_inv * S * J_inv;
                        v_HAC = J_inv * S_HAC * J_inv;
                        v_theta(i,j) = v(2,2); v_theta_HAC(i,j) = v_HAC(2,2);
                        v_all(((i-1)*(2)+1):(i*(2)),((j-1)*(2)+1):(j*(2))) = v;
                        v_all_HAC(((i-1)*(2)+1):(i*(2)),((j-1)*(2)+1):(j*(2))) = v_HAC;
                    end
                end
            end
        end
        function [T_out, p_out] = max_test_pval_fcn(obj, V, M, W, stat_in)
            k_delta = 1;
            %temp_v = chol(V);
            temp_v = sqrtm(V);
            Z = randn(obj.k_theta_n*(k_delta + 1),M);
            Z_sim = temp_v * Z;
            ind = (1+k_delta:1+k_delta:obj.k_theta_n*(k_delta + 1));
            Z_theta = Z_sim(ind,:);            
            W = repmat(W, 1, M);
            temp = abs(W .* Z_theta);
            for ind = 1:length(obj.k_theta_n_vec)
                statistic = stat_in{ind};
                k = obj.k_theta_n_vec(ind);
                temp2 = temp(1:k,:);
                T = max(temp2,[],1);
                T_out{ind} = T;
                p_out{ind} = obj.pval_fcn(T, M, statistic);
            end
       end
        function [T_out, p_out] = max_test_pval_bootstrap_fcn(obj, M, W_in, stat_in)
            X = ones(obj.n,1);
            beta_hat_null_imp = X \ obj.Y0; % Just OLS here;
            e_hat = obj.Y0 - beta_hat_null_imp' * X;
            z = randn(obj.n,M);    
            y_hat_null_imp_mat = beta_hat_null_imp' .* ones(obj.n,M) + (repmat(e_hat,1,M) .* z);
            theta_mat = zeros(obj.k_theta_n,M);
            for i = 1:obj.k_theta_n
                temp_x = [ones(obj.n-i,1) obj.Y0(1:end-i)];
                for m = 1:M
                    temp_y = y_hat_null_imp_mat(1+i:end, m);
                    b_hat = temp_x \ temp_y;
                    theta_mat(i,m) = b_hat(2);
                end
            end
            for ind = 1:length(obj.k_theta_n_vec)
                stat_in2 = stat_in{ind};
                k = obj.k_theta_n_vec(ind);
                T_out{ind} = zeros(M, length(stat_in2)); p_out{ind} = zeros(length(stat_in2),1);
                for q = 1:length(stat_in2)
                    W = W_in(:,q); statistic = stat_in2(q);
                    temp = abs(repmat(W,1,M) .* theta_mat);
                    temp2 = temp(1:k,:);
                    T = sqrt(obj.n) * max(temp2)';
                    T_out{ind}(:,q) = T;
                    p_out{ind}(q) = obj.pval_fcn(T, M, statistic);
                end
            end
        end
        
        % Ljung-Box Test
        function LBQ_vec = LBQ_fcn(obj, rho_vec)
            LBQ_vec = zeros(1,length(obj.k_theta_n_vec));
            rho_vec = rho_vec(1:obj.k_theta_n);
            for h = 1:obj.k_theta_n
                rho_vec(h) = ((obj.n-2)/(obj.n-h)) * (obj.n * (rho_vec(h)^2) - 1);
            end
            for ind = 1:length(obj.k_theta_n_vec)
                k = obj.k_theta_n_vec(ind);
                LBQ_vec(ind) = ((2 * k)^(-1/2)) * sum(rho_vec(1:k));
            end
        end
        
        % AP Sup-LM Test
        function sup_LM = sup_LM_fcn(obj, rho_vec)
            len_r = length(rho_vec);
            len_l = length(obj.lambda_grid);
            index = (1:len_r)';
            index = repmat(index, 1, len_l);
            lambdas = repmat(obj.lambda_grid, len_r, 1);
            lambdas = lambdas .^(index-1);
            rhos = repmat(rho_vec, 1, len_l);
            lambda_rho = lambdas .* rhos;
            LM = obj.n * (1-obj.lambda_grid.^2) .* sum(lambda_rho,1).^2;
            sup_LM = max(LM);
        end
        
        % CvM Test
        function CvM = CvM_fcn(obj, gamma_vec)
            len = length(obj.pi_grid);
            h_vec = (1:(obj.n-1))';
            psi_mat = repmat((1 ./ (h_vec .* pi)),1,len) .* sin(h_vec * obj.pi_grid);
            S = sqrt(obj.n) * gamma_vec' * psi_mat;
            CvM = sum(S.^2) / floor(len/pi);
        end
        
        % Max Correlation Test
        function T_max_corr_vec = max_corr_fcn(obj, rho_vec)
            T_max_corr_vec = zeros(1,length(obj.k_theta_n_vec));
            for ind = 1:length(obj.k_theta_n_vec)
                k = obj.k_theta_n_vec(ind);
                temp_rho_vec = rho_vec(1:k);
                T_max_corr_vec(ind) = sqrt(obj.n) * max(abs(temp_rho_vec));
            end
        end
        
        % General Functions
        function p = pval_fcn(obj, T, M, statistic)
            p = (1/M) * sum((T>statistic));
        end
        function gamma = covar_fcn(obj, h, w)
            h = abs(h);
            y = obj.Y0 - sum(obj.Y0) / obj.n;
            E = y(1:(end - h)) .* y((1 + h):end);
            if nargin <= 2
                gamma = sum(E) / obj.n;
            elseif nargin == 3
                gamma = sum( w((1 + h):end) .* (E - (sum(E) / obj.n)) ) / obj.n;
            end
        end
        function rho = corr_fcn(obj, h, w)
            if nargin <= 2
                gamma_h = obj.covar_fcn(h);    
            elseif nargin == 3
                gamma_h = obj.covar_fcn(h, w);
            end
            gamma_0 = obj.covar_fcn(0);
            rho = gamma_h / gamma_0;
        end
        function [gamma_vec, rho_vec] = fcn_rho_vec(obj, w)
            rho_vec = zeros((obj.n-1),1);
            if nargin <= 1
                for h = 1:(obj.n - 1)
                    rho_vec(h) = obj.corr_fcn(h);
                end
            elseif nargin == 2
                for h = 1:(obj.n - 1)
                    rho_vec(h) = obj.corr_fcn(h,w);
                end
            end
            gamma_vec = zeros((obj.n-1),1);
            if nargin <= 1
                for h = 1:(obj.n-1)  
                    gamma_vec(h) = obj.covar_fcn(h);
                end
            elseif nargin == 2
                for h = 1:(obj.n-1)  
                    gamma_vec(h) = obj.covar_fcn(h,w);
                end
            end            
        end
        % Dependent Wild Bootstrap
        function w = DWB_w_fcn(obj)
            num_windows = ceil((obj.n) / obj.bn);
            w1 = randn(1,num_windows);
            w1 = repmat(w1,obj.bn,1);
            w1 = w1(:);
            w = w1(1:obj.n);
        end
        function wM = DWB_wM_fcn(obj) % not used
            wM = zeros(obj.n, obj.M);
            for m = 1:obj.M
                wM(:,m) = obj.DWB_w_fcn;
            end
        end
    end
end
