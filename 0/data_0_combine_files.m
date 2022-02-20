% Compile data

clear all;
clc;

% You must set the correct working directory so the data has access to its
% class file
Test_number = 1; 
n_vec = [100 250 500 1000];

temp = pwd;
if Test_number ~= str2num(temp(end)) 
    fprintf('\n WARNING: Incorrect Working Directory - Will add to path.\n \n');
end
test_num_dir = sprintf('../%d', Test_number);
addpath(test_num_dir);


% data_main_dir = sprintf('G:/Simulation_data/max_test_many_zeros/kd_data/%d/data_n%d', Test_number, n);
data_main_dir0 = '../data';
data_out_dir0 = '../data_combined';

for n = n_vec
    data_dir0 = sprintf('%d/data_n%d', Test_number, n);
    data_main_dir = sprintf('%s/%s', data_main_dir0, data_dir0);

    if exist(data_main_dir, 'dir') == 7
        data_out_dir = sprintf('%s/%s', data_out_dir0, data_dir0);
        dir_list = dir(data_main_dir);
        N_dirs = length(dir_list);

        for d = 1:N_dirs
            if dir_list(d).isdir == 1 
                data_sub_dir = dir_list(d).name;
                if (strcmp(data_sub_dir,'.') == 0 && strcmp(data_sub_dir,'..') == 0)
                    data_dir = sprintf('%s/%s', data_main_dir, data_sub_dir);
                    outputdir = sprintf('%s/%s', data_out_dir, data_sub_dir);
                    fprintf('\n %d/%d:  %s \n', d, N_dirs, data_sub_dir);

                    file_list = dir(data_dir);

                    data_all = [];
                    N_files = length(file_list);
                    for file = 1:N_files
                        if file_list(file).isdir == 0
                            temp = file_list(file).name;
                            %if str2num(temp(12:end-4)) > 10
                            if temp(1:4) ~= 'comb'
                                file_name = sprintf('%s/%s',data_dir,temp);
                                load(file_name);
                                fprintf('%s \n', temp);
                                data_all = [data_all data];
                                clear data time
                            end
                            %end
                        end
                    end

                    data = data_all;
                    clear data_all;

                    if exist(outputdir,'dir') == 0
                        mkdir(outputdir);
                    end;

                    outputname=sprintf('%s/combined_%s.mat', outputdir, data_sub_dir);
                    save(outputname, 'data');
                    fprintf('NEW FILE: %s \n', outputname);

                end
            end
        end
    end % dir exists check
end % n_vec

rmpath(test_num_dir);

