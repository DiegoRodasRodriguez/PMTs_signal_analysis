close all
clear all

disp('                                                                       ')
disp('           *---------------------------------------------------------* ')
disp('           *          Macro for reading several ASCII files from     * ')
disp('           *                     MADX and convert to .mat            * ')
disp('           *---------------------------------------------------------* ')
disp('                                                                       ')


Nfiles    = 10;
Directory = 'E:\MATLAB7\work\Angela\';
File      = 'test_header';
Ext       = 'txt';
Ncol      = 10;

Dir_var    = 'Directory';
File_var   = 'File';

for i=1:Nfiles    
    Ext_    = [Ext num2str(i)];
    Ext_var = 'Ext_';
    
    eval(['readrun_A(' Dir_var ',' File_var ',' Ext_var ',' num2str(Ncol) ');']);    
end
   