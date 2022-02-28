close all; clear all; clc;

% types = {'Ca cha'; 'PSD95'; 'syngap'; 'homer'; 'Rim'; 'basoon'; 'shank'};
% layers = {'L1'; 'L23'; 'L4'};
% AB = {'A'; 'B'};

filename = 'result.xlsx';
sheetname_set = {'n=7'};

% iterate over sheets
for sheet_idx = 1:size(sheetname_set,1)
    sheetname = sheetname_set{sheet_idx};
    
    result = readtable(filename, 'Sheet', sheetname);
    
    figure; 
    bar(result{:,1}, result{:,2});
    errorbar(result{:,2}, result{:,4});
    
end



