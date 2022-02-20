
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



data_main_dir0 = '../data_combined';

for n = n_vec
    data_dir0 = sprintf('%d/data_n%d', Test_number, n);
    data_main_dir = sprintf('%s/%s', data_main_dir0, data_dir0);

    if exist(data_main_dir, 'dir') == 7
        dir_list = dir(data_main_dir);
        N_dirs = length(dir_list);

        for d = 1:N_dirs
            if dir_list(d).isdir == 1 
                data_sub_dir = dir_list(d).name;
                if (strcmp(data_sub_dir,'.') == 0 && strcmp(data_sub_dir,'..') == 0)
                    data_dir = sprintf('%s/%s', data_main_dir, data_sub_dir);
                    fprintf('\n %d/%d:  %s \n', d, N_dirs, data_sub_dir);

                    file_list = dir(data_dir);
                    N_files = length(file_list);
                    for file = 1:N_files
                        if file_list(file).isdir == 0
                            temp = file_list(file).name;
                            if temp(1:4) == 'comb'  %only use the combined data file%
                                file_name = sprintf('%s/%s',data_dir,temp);
                                load(file_name);
                                fprintf('%s \n', temp);
                                %
                                data_1_1_post_process;
                                %
                            end
                        end
                    end
                    clearvars -except dir_list N_dirs data_main_dir d test_num_dir Test_number data_main_dir0;
                end
            end
        end
    end % dir exists check
end % n_vec
rmpath(test_num_dir);
