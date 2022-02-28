function [volume, SNR, protein, layerOrder, ABOrder, prePost] = readResultForANOVA(...
    sheetname, types, layers, AB)

volume = [];
SNR = [];
protein = {};
layerOrder = {};
ABOrder = {};
prePost = {};

% read from directory
resultFile = dir('result/*.xlsx');

% iterate over result excel files
for file_idx = 1:size(resultFile,1)
    seqname = resultFile(file_idx).name;
    result = readtable(fullfile('result', seqname), 'Sheet', sheetname);
    
    % read one row in one excel file
    for data_idx = 1:size(result, 1) - 6 % exclude average value
        
        % read/write layer order
        for layer_idx = 1:size(layers,1)
            if contains(seqname, cell2mat(layers(layer_idx)))
                layerOrder = [layerOrder; layers(layer_idx); layers(layer_idx)];
                break;
            end
            if layer_idx == size(layers,1)
                disp("ERROR: TYPE DOES NOT EXISTS"); break;
            end
        end
                
        % read/write AB order
        for AB_idx = 1:size(AB,1)
            if contains(seqname(1), cell2mat(AB(AB_idx)))
                ABOrder = [ABOrder; AB(AB_idx); AB(AB_idx)];
                break;
            end
            if AB_idx == size(AB,1)
                disp("ERROR: TYPE DOES NOT EXISTS"); break;
            end
        end
        
        % read/write protein order
        if size(seqname) < 10
            len_dirname = size(seqname);
        else
            len_dirname = 10;
        end
        for type_idx = 1:size(types,1)
            if contains(seqname(1:len_dirname), cell2mat(types(type_idx)))
                protein = [protein; types(type_idx); types(type_idx)];
                break;
            end
            if type_idx == size(types,1)
                disp("ERROR: TYPE DOES NOT EXISTS"); break;
            end
        end
        
        
        % read/write pre/post order
        prePost = [prePost; {'pre'}; {'post'}];
        

        
        
        %% read/write volume, SNR values

        % volume, pre
        if iscell(result{data_idx, 2})
            volume = [volume; str2double(result{data_idx, 2}{:})]; % pre
        else
            volume = [volume; result{data_idx, 2}]; % pre
        end
        
        % volume, post
        if iscell(result{data_idx, 3})
            volume = [volume; str2double(result{data_idx, 3}{:})]; % post
        else
            volume = [volume; result{data_idx, 3}]; % post
        end
        
        % SNR, pre
        if iscell(result{data_idx, 4})
            SNR = [SNR; str2double(result{data_idx, 4}{:})]; % pre
        else
            SNR = [SNR; result{data_idx, 4}]; % pre
        end
        
        % SNR, post
        if iscell(result{data_idx, 5})
            SNR = [SNR; str2double(result{data_idx, 5}{:})]; % post
        else
            SNR = [SNR; result{data_idx, 5}]; % post
        end
    end
end
end