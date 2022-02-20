if exist('data', 'var') == 1

    J = length(data);
    rng_seed = [];
    for j = 1:J
        rng_seed = [rng_seed; data(j).rng_seed];
    end

    unique_rng = length(unique(rng_seed));
    if unique_rng ~= J
        fprintf('\n \n ERROR: Unique Seeds = %d, but J=%d \n \n ', unique_rng, J);
        pause;
    end

    %rej_table;
    
    clear rej_table temp_dr;
    J = length(data);
    temp_dr = zeros(size(data(1).dr,1),size(data(1).dr,2),J);
    for j = 1:J
        temp_dr(:,:,j) = data(j).dr;
    end
    rej_table = sum(temp_dr,3)/J;


    %% Make tex table
    outputdir = '../tables';
    test_number = data(1).test_number;  

    if exist(outputdir,'dir') == 0
        mkdir(outputdir);
    end
    table_name = 'RF';

    clear table1 rownames colnames Title;

    if test_number == 1
        rownames{1} = sprintf('Max Test');
        rownames{2} = sprintf('Max t-Test');
        rownames{3} = sprintf('Max t-Test (HAC)');
        rownames{4} = sprintf('Max Test (PBS)');
        rownames{5} = sprintf('Max t-Test (PBS)');
        rownames{6} = sprintf('Max t-Test (HAC) (PBS)');
        rownames{7} = sprintf('LBQ');
        rownames{8} = sprintf('sup LM');
        rownames{9} = sprintf('CvM');
        rownames{10} = sprintf('Max Corr');
    elseif test_number == 2
        rownames{1} = sprintf('Max Test');
        rownames{2} = sprintf('Max t-Test');
        rownames{3} = sprintf('Max Test (BS)');
        rownames{4} = sprintf('Max t-Test (BS)');
        rownames{5} = sprintf('Wald');
        rownames{6} = sprintf('Wald (BS)');
    elseif test_number == 3
        rownames{1} = sprintf('Max Test');
        rownames{2} = sprintf('Max t-Test');
        rownames{3} = sprintf('LBQ');
        rownames{4} = sprintf('Max Corr');
        rownames{5} = sprintf('AST');
        rownames{6} = sprintf('FZ Wald');
        rownames{7} = sprintf('FZ Wald (NS)');
    end

    alpha_levels = data(1).alpha_levels;
    for a = 1:length(alpha_levels)
        colnames{a} = sprintf('$(\\alpha=%.2f)$', alpha_levels(a));
    end

    switch test_number
        case 1
            dgp_type = data(1).dgp_type; dgp_type_string = data(1).dgp_type_string;
            % hypothesis_type = data(1).hypothesis_type; hypothesis_type_string = data(1).hypothesis_type_string;
            n = data(1).n;
            % k_theta = data(1).k_theta;
            % k_delta = data(1).k_delta;
            k_theta_n = data(1).k_theta_n;
            Caption=sprintf('Rejection Frequencies: %s, $n=%d$, $k_{\\theta,n}=%d$', dgp_type_string, n, k_theta_n);
            Title{1}=sprintf('Rejection Frequencies, J=%d', J);
            Title{2}=sprintf('Experiment: %d, DGP: %s', test_number, dgp_type_string);
            Title{3}=sprintf('$n=%d$, $k_{\\theta,n}=%d$', n, k_theta_n);
            outputname = sprintf('./%s/%s_%d_dgp%d_n%d_ktn%d', outputdir, table_name, test_number, dgp_type, n, k_theta_n);
        case 2
            dgp_type = data(1).dgp_type; dgp_type_string = data(1).dgp_type_string;
            hypothesis_type = data(1).hypothesis_type; hypothesis_type_string = data(1).hypothesis_type_string;
            n = data(1).n;
            k_theta = data(1).k_theta;
            k_delta = data(1).k_delta;
            k_theta_n = data(1).k_theta_n;
            Caption=sprintf('Rejection Frequencies: %s, %s, $n=%d$, $k_{\\delta}=%d$, $k_{\\theta}=%d$, $k_{\\theta,n}=%d$', dgp_type_string, hypothesis_type_string, n, k_delta, k_theta, k_theta_n);
            Title{1}=sprintf('Rejection Frequencies, J=%d', J);
            Title{2}=sprintf('Experiment: %d, DGP: %s, Hyp: %s', test_number, dgp_type_string, hypothesis_type_string);
            Title{3}=sprintf('$n=%d$, $k_{\\delta}=%d$, $k_{\\theta}=%d$, $k_{\\theta,n}=%d$', n, k_delta, k_theta, k_theta_n);
            outputname = sprintf('./%s/%s_%d_dgp%d_hyp%d_n%d_kd%d_kt%d_ktn%d', outputdir, table_name, test_number, dgp_type, hypothesis_type, n, k_delta, k_theta, k_theta_n);
        case 3
            dgp_type = data(1).dgp_type; dgp_type_string = data(1).dgp_type_string;
            % hypothesis_type = data(1).hypothesis_type; hypothesis_type_string = data(1).hypothesis_type_string;
            n = data(1).n;
            % k_theta = data(1).k_theta;
            % k_delta = data(1).k_delta;
            k_theta_n = data(1).k_theta_n;
            Caption=sprintf('Rejection Frequencies: %s, $n=%d$, $k_{\\theta,n}=%d$', dgp_type_string, n, k_theta_n);
            Title{1}=sprintf('Rejection Frequencies, J=%d', J);
            Title{2}=sprintf('Experiment: %d, DGP: %s', test_number, dgp_type_string);
            Title{3}=sprintf('$n=%d$, $k_{\\theta,n}=%d$', n, k_theta_n);
            outputname = sprintf('./%s/%s_%d_dgp%d_n%d_ktn%d', outputdir, table_name, test_number, dgp_type, n, k_theta_n);
    end

    table1 = rej_table;

    tabletotex(table1, rownames, colnames, outputname, Title);
 
end

