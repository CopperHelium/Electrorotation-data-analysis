%----------------------------------------------------------------
%Zou Ying - v1.0
%
%----------------------------------------------------------------
%Navish Wadhwa - v1.5
%This area contians are the parameter which may be changed every times
path='E:\work\08-04\';

%----------------------------------------------------------------
Files=dir(strcat(path,'*.tdms'));
lengthFiles=length(Files);
for i=1:lengthFiles
   % fopen(strcat('E:\work\',Files(i).name),'rt');
    simpleConvertTDMS(strcat(path,Files(i).name));
end

fprintf('Completed''\n');
clear;

