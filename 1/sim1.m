function sim1(dgp_type, n, J, seed_type)

if nargin == 0
    clear all; clc;
    J = 2;            %number of Simulations
    dgp_type = 1;
    n = 99;
    seed_type = 1;
end

test_number = 1;

% dgp_type_vec = 1:7;
% n_vec = [100, 250, 500, 1000];
k_theta_n_vec = [5, floor(3*n^(2/5)), floor(7*n^(2/5))];
if (dgp_type == 7 && n == 100)
    k_theta_n_vec = [k_theta_n_vec, 50];
end

warning('off','all');

if seed_type == 0  % all at once, random
    seed_input = 'shuffle';
    rng('shuffle');
    temp = rng;
    sim_number = temp.Seed;
    clear temp;
elseif seed_type == 1  % all at once, reproducible
    seed_input0 = 0;
    sim_number = seed_input0;
elseif seed_type == 2  % longleaf array, reproducible
    seed_input0 = str2num(getenv('SLURM_ARRAY_TASK_ID'));
    seed_input0 = seed_input0 - 1;
    sim_number = seed_input0;
end

data0 = [];

tic;
for j = 1:J
    if (seed_type == 1 || seed_type == 2)
        seed_input = seed_input0 * 10^ceil(log(J)/log(10)) + j;
    end
    data0 = [data0 class_tests_1(dgp_type, n, seed_input)];
end
time.dgp = toc / 60


    tic;
    for j = 1:J
        data0(j) = data0(j).run_all_tests_fcn(k_theta_n_vec);
    end
    time.est = toc / 60

    %clean up
    for j = 1:J
        data0(j) = data0(j).clean_up();
    end

for k_theta_n = k_theta_n_vec
    data = data0;
    for j = 1:J
        data(j) = data(j).fcn_final_output(k_theta_n);
    end

    output_main_dir = sprintf('./data_n%d', n);

    outputdir = sprintf('%s/output_%d_dgp%d_n%d_ktn%d', output_main_dir, test_number, dgp_type, n, k_theta_n);
    if exist(outputdir,'dir') == 0
        mkdir(outputdir);
    end
    outputname=sprintf('%s/data_%d_dgp%d_n%d_ktn%d_J%d_%d.mat', outputdir, test_number, dgp_type, n, k_theta_n, J, sim_number);

    sim_number0 = sim_number;
    while exist(outputname,'file') == 2
        display('Error: Name Taken; Iterating Name Forward');
        sim_number0 = sim_number0 + 1;
        outputname=sprintf('%s/data_%d_dgp%d_n%d_ktn%d_J%d_%d.mat', outputdir, test_number, dgp_type, n, k_theta_n, J, sim_number0);
    end
    display(outputname);
    save(outputname, 'data', 'time');

    clear data;

end % k_theta_n loop

clear data0 time;

% end % dgp_type loop

end % function

