classdef class_dgp_1
    properties (SetAccess=public, GetAccess=public)
	 %dgp variables
     test_number
     dgp_type
     dgp_type_string
     n
     init_n
	 Y0
     e
     seed_input
     rng_seed
     %estimation variables
    end
    methods
        function obj = class_dgp_1(dgp_type, n, seed_input)
            obj.test_number = 1;
            obj.dgp_type = dgp_type;
            obj.n = n;
            obj.seed_input = seed_input;
            rng(seed_input);
            temp = rng;
            obj.rng_seed = temp.Seed;
            clear temp;
            % Now generate the data
	        obj.init_n = n + n;
            obj.e = randn(obj.init_n,1);
            switch dgp_type
                case 1
                    obj.Y0 = obj.dgp_iid;
                    obj.dgp_type_string = 'iid';
                case 2
                    obj.Y0 = obj.dgp_GARCH;
                    obj.dgp_type_string = 'GARCH';
                case 3
                    obj.Y0 = obj.dgp_bilinear;
                    obj.dgp_type_string = 'Bilinear';
                case 4
                    obj.Y0 = obj.dgp_AR1;
                    obj.dgp_type_string = 'AR(1)';
                case 5
                    obj.Y0 = obj.dgp_MA1;
                    obj.dgp_type_string = 'MA(1)';
                case 6
                    obj.Y0 = obj.dgp_MA20;
                    obj.dgp_type_string = 'MA(20)';
                case 7
                    obj.Y0 = obj.dgp_MA50;
                    obj.dgp_type_string = 'MA(50)';
            end
            obj.Y0 = obj.Y0((obj.init_n - n + 1):end); % remove Y0 burn-in values
        end
        function y = dgp_iid(obj)
            y = obj.e;
        end
        function y = dgp_GARCH(obj,beta)
            beta = [.3, .6];
            sigma2 = ones(obj.init_n,1);
            y = zeros(obj.init_n,1);
            y(1) = sqrt(sigma2(1)) * obj.e(1);
            for t = 2:obj.init_n
                sigma2(t) = 1 + beta(1) * y(t-1)^2 + beta(2) * sigma2(t-1);
                y(t) = sqrt(sigma2(t)) * obj.e(t);
            end
        end
        function y = dgp_bilinear(obj,beta)
            beta = 1;
            y = zeros(obj.init_n,1);
            y(1) = obj.e(1);
            for t = 2:obj.init_n
                y(t) = beta * y(t-1) * obj.e(t-1) + obj.e(t);
            end
        end
        function y = dgp_AR1(obj, beta)
            beta = .5;
            y = zeros(obj.init_n,1);
            for t = 2:obj.init_n
                y(t) = beta * y(t-1) + obj.e(t);
            end
        end
        function y = dgp_MA1(obj,beta)
            beta = 1;
            y = zeros(obj.init_n,1);
            for t = 2:obj.init_n
                y(t) = beta * obj.e(t-1) + obj.e(t);
            end
        end
        function y = dgp_MA20(obj,beta)
            beta = 1;
            y = zeros(obj.init_n,1);
            for t = 21:obj.init_n
                y(t) = beta * obj.e(t-20) + obj.e(t);
            end
        end
        function y = dgp_MA50(obj,beta)
            beta = 1;
            y = zeros(obj.init_n,1);
            for t = 51:obj.init_n
                y(t) = beta * obj.e(t-50) + obj.e(t);
            end
        end
    end
end
