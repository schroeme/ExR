close all; clear all; clc;

types = {'Ca cha'; 'PSD95'; 'syngap'; 'homer'; 'Rim'; 'basoon'; 'shank'};
layers = {'L1'; 'L23'; 'L4'};
AB = {'A'; 'B'};

filename = 'result_calculated_Feb_26.xlsx';
sheetname_set = {'n=7'};

repeatN = size(AB,1); % AB

% iterate over sheets
for sheet_idx = 1:size(sheetname_set,1)
    
    sheetname = sheetname_set{sheet_idx};
    result = readtable(filename, 'Sheet', sheetname);

    [volume, SNR, protein, layerOrder, prePost] = readResultForANOVA(...
        result, repeatN, types, layers);
    
    [p_vol, tbl_vol, stats_vol] = anovan(volume, {protein, layerOrder, prePost}, ...
        'varnames',{'protein','layer','prePost'});
%         'model','interaction','varnames',{'protein','layer','prePost'});
    [p_SNR, tbl_SNR, stats_SNR] = anovan(SNR, {protein, layerOrder, prePost}, ...
        'varnames',{'protein','layer','prePost'});
%         'model','interaction','varnames',{'protein','layer','prePost'});

    
end







