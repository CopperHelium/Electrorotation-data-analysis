%----------------------------------------------------------------
%Zou Ying - v1.0
%
%----------------------------------------------------------------
%Navish Wadhwa - v1.5
% - Made user select folder. 
% - updated comments. 
% - general cleanup.

%Ask the user to select folder. use current working directory as a starting
%point. 

start_path=pwd;
folder_name = uigetdir(start_path,'Select folder containing tdms files to be converted to mat:');
cd(folder_name)
%----------------------------------------------------------------

% Collect a list of all the files in the current directory
Files=dir('*.tdms');

% Measure how long the list of files is
lengthFiles=length(Files);

%Use simpleConvertTDMS function to convert all files into mat files one
%after the other. All mat files will be saves in the same directory as the
%tdms files
for i=1:lengthFiles
       simpleConvertTDMS(Files(i).name);
end

% Display completion message when done
disp('All done.')


