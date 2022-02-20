% data 2: output latex file with includes for all tables and figures.
clear all; clc;

main_dir = '..';
table_dir = 'tables';
fig_dir = 'figs';
fig_scale = '.5';
output_name = 'max_test_many_zeros';
outputdir = '../tables_combined';
if exist(outputdir,'dir') == 0
    mkdir(outputdir);
end

%% Tables
data_dir = sprintf('%s/%s', main_dir, table_dir);
if exist(data_dir,'dir') == 7
    file_list = dir(data_dir);
    fname = sprintf('%s/tables_%s.tex', outputdir, output_name);
    FID = fopen(fname, 'w');
%     fprintf(FID, '\\documentclass[11pt]{article} \n');
%     fprintf(FID, '\\usepackage{setspace,graphicx,caption,subcaption,float,relsize} \n');
    fprintf(FID, '\\input{./sections/header_article.tex} \n');
    fprintf(FID, '\\input{./sections/header_commands.tex} \n \n');
    fprintf(FID, '\\begin{document} \n \n');
    fprintf(FID, '\\title{{\\large Supplemental Appendix - Monte Carlo Simulations for \\medskip \\\\ Testing Many Zero Restrictions \\medskip \\\\ Where a Subset May Lie On the Boundary }} \n \n');
    fprintf(FID, '\\author{ \\ \\ Jonathan B. Hill\\thanks{Department of Economics, University of North Carolina, Chapel Hill. E-mail: \\texttt{jbhill@email.unc.edu}} : University of North Carolina \\medskip \\\\ Jay Dennis\\thanks{E-mail: \\texttt{jay.dennis@unc.edu}} : University of North Carolina} \n \n');
    fprintf(FID, '\\date{{\\large This draft:} \\today } \n \n');
    fprintf(FID, '\\maketitle \n \n');
    fprintf(FID, '\\newpage \n \n');
    fprintf(FID, '\\tableofcontents \n \n');
    %
    % Insert code here
    table_name = 'RF';
    n_vec = [100 250 500 1000];
    test_number_vec = 1:3;
    for test_number = test_number_vec
        fprintf(FID, '\\newpage \n \n');
        fprintf(FID, '\\section{Experiment %d} \n \n', test_number);                    
            fname_obs = sprintf('./Observations_%d.txt', test_number);
            if exist(fname_obs,'file') == 2
                temp_file_ID = fopen(fname_obs, 'r+');
                while 1
                    temp_line = fgets(temp_file_ID);
                    if ~ischar(temp_line) 
                        fprintf(FID, '\n \n');
                        break
                    end
                    fprintf(FID, '%s', temp_line);
                    %fprintf('%s', temp_line);
                end                
                fclose(temp_file_ID);
            end
        fprintf(FID, '\\newpage \n \n');
        if test_number == 1
            dgp_type_vec = 1:7;
            hypothesis_type_vec = 1;
        elseif test_number == 2
            dgp_type_vec = 1:3;
            hypothesis_type_vec = 1:5;
        elseif test_number == 3
            dgp_type_vec = 1:7;
            hypothesis_type_vec = 1;
        end
        for dgp_type = dgp_type_vec
            for n = n_vec
                if test_number == 1
                    k_theta_n_vec = unique([5, floor(3*n^(2/5)), floor(7*n^(2/5)), 50]);
                    k_theta_vec = 1;
                    k_delta_vec = 1;
                elseif test_number == 2
                    k_theta_vec = [10, 40, floor(n/2)];
                    k_delta_vec = [0, 20];
                    k_theta_n_vec = unique([1, 10, 20, floor(5*n^(1/2 - .0001))]);
                elseif test_number == 3
                    k_theta_n_vec = unique([1, 10, 20, floor(25*n^(.16)), floor(2*(log(n))^2)]);
                    k_theta_vec = 1;
                    k_delta_vec = 1;
                end
                for hypothesis_type = hypothesis_type_vec
                    if test_number ~= 2
                        fprintf(FID, '\\subsection{Experiment %d, DGP %d, n = %d} \n \n', test_number, dgp_type, n);
                    elseif test_number == 2
                        fprintf(FID, '\\subsection{Experiment %d, DGP %d, Hypothesis %d, n = %d} \n \n', test_number, dgp_type, hypothesis_type, n);
                    end
                    for k_theta = k_theta_vec
                        for k_delta = k_delta_vec
                            for k_theta_n = k_theta_n_vec
                                if test_number == 1
                                    temp = sprintf('%s/%s_%d_dgp%d_n%d_ktn%d', data_dir, table_name, test_number, dgp_type, n, k_theta_n);
                                elseif test_number == 2
                                    temp = sprintf('%s/%s_%d_dgp%d_hyp%d_n%d_kd%d_kt%d_ktn%d', data_dir, table_name, test_number, dgp_type, hypothesis_type, n, k_delta, k_theta, k_theta_n);
                                elseif test_number == 3
                                    temp = sprintf('%s/%s_%d_dgp%d_n%d_ktn%d', data_dir, table_name, test_number, dgp_type, n, k_theta_n);
                                end
                                fprintf('%s \n', temp);
                                %
                                file_in = sprintf('%s.tex', temp);
                                if exist(file_in,'file') == 2
                                    temp_file_ID = fopen(file_in, 'r+');
                                    while 1
                                        temp_line = fgets(temp_file_ID);
                                        if ~ischar(temp_line) 
                                            fprintf(FID, '\n \n');
                                            break
                                        end
                                        fprintf(FID, '%s', temp_line);
                                        %fprintf('%s', temp_line);
                                    end
                                    fclose(temp_file_ID);
                                end
                                %               
                            end
                        end
                    end
                end
            end
        end
    end
    %
    fprintf(FID, '\\end{document} \n');
    fclose(FID);

end


