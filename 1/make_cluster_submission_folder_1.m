clear all;
clc;

cluster_job_type = 'array';

sim_file_name = 'sim1';
test_number = 1;
dgp_type_vec = (1:7);
T_vec = [100 250 500 1000];

starting_sim_number = 5001; % be careful here. 
total_num_sims = 5000;
num_sims_per_array_job = 100;
num_array_jobs = ceil(total_num_sims / num_sims_per_array_job);
array_start = floor(starting_sim_number / num_sims_per_array_job) + 1;
array_len{100} = num_array_jobs; array_len{250} = num_array_jobs; array_len{500} = num_array_jobs; array_len{1000} = num_array_jobs;

array_batch_size = 3;
J_vec{100} = total_num_sims; J_vec{250} = total_num_sims; J_vec{500} = total_num_sims; J_vec{1000} = total_num_sims;

infile = fullfile(sprintf('%s.m', sim_file_name));

cluster_dir = sprintf('./cluster_sub');
if exist(cluster_dir,'dir') == 0
    mkdir(cluster_dir);
end

for T = T_vec
    if strcmp(cluster_job_type, 'array') == 1
        J = floor(total_num_sims / array_len{T});
        array_end = array_start + array_len{T} - 1;
    elseif strcmp(cluster_job_type, 'all') == 1
        J = J_vec{T};
    end
    outputname_ll = sprintf('ll_%d_T%d', test_number, T);
    fname = sprintf('%s/%s.txt', cluster_dir, outputname_ll);
    FID = fopen(fname, 'w');
    %allocate 5 min/sim (for longleaf) * J * T/100
    wall_time = 5*J*T/100;
    for dgp_type = dgp_type_vec
        job_name = sprintf('MT%dd%dT%dJ%d', test_number, dgp_type, T, J);
        if strcmp(cluster_job_type, 'array') == 1
            m_file_name = sprintf('%s(%d,%d,%d,2)', sim_file_name, dgp_type, T, J);
            fprintf(FID, 'sbatch --job-name=%s -o %s.%%j.txt -t %d --array=[%d-%d%%%d] --wrap=\"matlab -nodisplay -nojvm -nodesktop -nosplash -singleCompThread -r ''%s; quit;''\"; \n \n', job_name, job_name, wall_time, array_start, array_end, array_batch_size, m_file_name);
        elseif strcmp(cluster_job_type, 'all') == 1
            m_file_name = sprintf('%s(%d,%d,%d,1)', sim_file_name, dgp_type, T, J);
            fprintf(FID, 'sbatch --job-name=%s -o %s.%%j.txt -t %d --wrap=\"matlab -nodisplay -nojvm -nodesktop -nosplash -singleCompThread -r ''%s; quit;''\"; \n', job_name, job_name, wall_time, m_file_name);
        end
        fprintf(FID, '\n');
    end
    fclose(FID);
end

outfile = fullfile(cluster_dir, sprintf('%s.m', sim_file_name));
copyfile(infile, outfile);

tempfile = sprintf('class_dgp_%d.m', test_number);
infile = fullfile(tempfile);
outfile = fullfile(cluster_dir, tempfile);
copyfile(infile, outfile);

tempfile = sprintf('class_tests_%d.m', test_number);
infile = fullfile(tempfile);
outfile = fullfile(cluster_dir, tempfile);
copyfile(infile, outfile);


