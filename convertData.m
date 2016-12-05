%----------------------------------------------------------------
%Zou Ying - v1.0
%
%----------------------------------------------------------------
%Navish Wadhwa - v1.5

%This area contians the path to the folder containing the tdms files to be
%converted. This should be changed everytime new files are to be converted.

path='E:\work\08-04\';

%----------------------------------------------------------------

% Collect a list of all the files in the current directory
Files=dir(strcat(path,'*.tdms'));

% Measure how long the list of files is
lengthFiles=length(Files);

%Use simpleConvertTDMS function to convert all files into mat files one
%after the other. All mat files will be saves in the same directory as the
%tdms files
for i=1:lengthFiles
       simpleConvertTDMS(strcat(path,Files(i).name));
end

% Display completion message when done
fprintf('Completed''\n');
clear;

