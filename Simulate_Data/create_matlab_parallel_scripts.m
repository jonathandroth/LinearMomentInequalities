

%% Automatically generate the Odyssey run scripts for the specified Matlab files using a Template
run_file_attributes = dir('run_lp*.m');


for filenum = 1:size(run_file_attributes)
    
  run_file_name = run_file_attributes(filenum).name;   
  run_file_name_noextension = strrep( run_file_name, '.m', '')
  
  
  fin = fopen('Matlab_parallel_script_template.sh');
  fout = fopen(strcat('Matlab_parallel_script_',run_file_name_noextension, '.sh'), 'w');
  
  while ~feof(fin)
      s = fgetl(fin);
      s = strrep(s, 'FILENAME', run_file_name_noextension);
      
      %Add more memory if 9thetacs
      if( ~isempty(regexp(run_file_name, '9thetac')) )
                s = strrep(s, '80000', '160000');
      end
      
      fprintf(fout,'%s\n',s);
      %disp(s)
  end
  
  fclose(fin);
  fclose(fout);

end