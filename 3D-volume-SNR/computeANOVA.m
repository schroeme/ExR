close all; clear all; clc;
%delete(findall(0));
 
proteinTypes = {'Ca cha'; 'PSD95'; 'syngap'; 'homer'; 'Rim'; 'basoon'; 'shank'};
layers = {'L1'; 'L23'; 'L4'};
AB = {'A'; 'B'};

sheetname_set = {'n=7'};

% iterate over sheets
for sheet_idx = 1:size(sheetname_set,1)
    sheetname = sheetname_set{sheet_idx};

    [volume, SNR, protein, layerOrder, ABOrder, prePost] = readResultForANOVA(...
        sheetname, proteinTypes, layers, AB);
    
    [p_vol, tbl_vol, stats_vol] = anovan(volume, {protein, layerOrder, prePost}, ...
        'varnames',{'protein', 'layer', 'prePost'});
%         'model','interaction','varnames',{'protein','layer','prePost'});

    [p_SNR, tbl_SNR, stats_SNR] = anovan(SNR, {protein, layerOrder, prePost}, ...
        'varnames',{'protein', 'layer', 'prePost'});
%         'model','interaction','varnames',{'protein','layer','prePost'});

    
end







