function F = dynare_cleanup( )

try
mod_files = dir('*.mod');
mod_files = cellstr(char({mod_files.name}));
num_mod_files = length( mod_files );

delete *_IRF_*
delete *_dynamic.m
delete *_static.m
delete *_steadystate2.m
delete *_results.mat
delete *.eps
delete *.pdf
delete *_set_auxiliary_variables.m
delete *.jnl
delete *.log
delete *Temp*
for i = 1:num_mod_files
    file_name = char(mod_files(i,:));
    if exist([file_name(1:end-4),'.m'], 'file')==2
        delete([file_name(1:end-4),'.m']);
    end
    if exist(file_name(1:end-4), 'dir')==7
        rmdir(file_name(1:end-4),'s');
    end
end
F = 1;
catch
F = 0;
end