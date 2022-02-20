classdef class_dgp_2
    properties (SetAccess=public, GetAccess=public)
	 %dgp variables
     test_number
     dgp_type
     dgp_type_string
     hypothesis_type
     hypothesis_type_string
     n
     init_n
	 Y0
     e
     k_delta
     k_theta
     x_delta
     x_theta
     delta_0
     theta_0
     seed_input
     rng_seed
     %estimation variables
    end
    methods
        function obj = class_dgp_2(dgp_type, hypothesis_type, n, seed_input, k_theta, k_delta)
            obj.test_number = 2;
            obj.dgp_type = dgp_type;
            obj.hypothesis_type = hypothesis_type;
            obj.n = n;
            obj.k_delta = k_delta;
            obj.k_theta = k_theta;
            obj.seed_input = seed_input;
            rng(seed_input);
            temp = rng;
            obj.rng_seed = temp.Seed;
            clear temp;
            % Now generate the data
	        obj.init_n = obj.n + 0;
            obj.e = randn(1,obj.init_n);
            switch dgp_type
                case 1
                    obj.dgp_type_string = 'Case 1: serial and multial iid';
                    obj.x_delta = randn(obj.k_delta,obj.n);
                    obj.x_theta = randn(obj.k_theta,obj.n);
                case 2 
                    obj.dgp_type_string = 'Case 2: blockwise independence';
                    w_delta = randn(obj.k_delta,obj.n);
                    w_theta = randn(obj.k_theta,obj.n);
                    rho = .5;
                    num = obj.k_delta;
                    Sig = eye(num) + rho .* (ones(num,num) - eye(num)); 
                    Sig = chol(Sig);
                    obj.x_delta = Sig * w_delta;
                    num = obj.k_theta;
                    Sig = eye(num) + rho .* (ones(num,num) - eye(num)); 
                    Sig = chol(Sig);
                    obj.x_theta = Sig * w_theta;
                case 3 
                    obj.dgp_type_string = 'Case 3: with-in and cross-block dependence';
                    num = obj.k_delta + obj.k_theta;
                    temp = randn(num,obj.n);
                    rho = .5;
                    Sig = eye(num) + rho .* (ones(num,num) - eye(num));
                    Sig = chol(Sig);
                    X = Sig * temp;
                    obj.x_delta = X(1:obj.k_delta,:);
                    obj.x_theta = X(obj.k_delta+1:end,:);
            end
            switch hypothesis_type
                case 1
                    obj.hypothesis_type_string = 'Null';
                    obj.theta_0 = zeros(obj.k_theta,1);
                case 2
                    obj.hypothesis_type_string = 'Alternative 1';
                    obj.theta_0 = zeros(obj.k_theta,1);
                    obj.theta_0(1) = 1/2;
                case 3
                    obj.hypothesis_type_string = 'Alternative 2';
                    obj.theta_0 = zeros(obj.k_theta,1);
                    for i = 1:10
                        obj.theta_0(i) = i/2;
                    end
                case 4
                    obj.hypothesis_type_string = 'Alternative 3';
                    obj.theta_0 = zeros(obj.k_theta,1);
                    for i = 1:obj.k_theta
                        obj.theta_0(i) = i/obj.k_theta;
                    end
                case 5
                    obj.hypothesis_type_string = 'Alternative 4';
                    obj.theta_0 = repmat(1/2, obj.k_theta, 1);
            end
            obj.delta_0 = ones(obj.k_delta,1);
            
            obj.Y0 = obj.delta_0' * obj.x_delta + obj.theta_0' * obj.x_theta + obj.e;
            
            obj.Y0 = obj.Y0((obj.init_n - obj.n + 1):end); % remove Y0 burn-in values
        end
    end
end
