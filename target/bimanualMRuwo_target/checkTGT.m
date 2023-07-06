clear;
close all;
clc;

addpath(genpath('/Users/aghavampour/Documents/MATLAB/dataframe-2016.1'),'-begin');

subjName = 'suwo06';

scan_file_name = strcat(subjName,'_scan_');

% going through runs:
for i = 1:10
    fname = strcat(scan_file_name,num2str(i),'.tgt');
    
    % loading .tgt file:
    tgt = datload(fname);

    % getting start times:
    start_time = tgt.startTime;

    % calculating the distances of strat times:
    diff_st = diff(start_time);

    % making sure of the values:
    [unique_st, ia, ic] = unique(diff_st);
    unique_st

    % counting the start times:
    for j = 1:length(unique_st)
        fprintf('diff = %d , occurance = %d\n',unique_st(j) , sum(ic == j));
    end

end

