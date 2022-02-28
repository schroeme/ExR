function [volume, SNR, protein, layerOrder, prePost] = readResultForANOVA(result, ...
    repeatN, types, layers)


volume = zeros(repeatN * size(result,1), 1);
SNR = zeros(repeatN * size(result,1), 1);
protein = cell(repeatN * size(result,1), 1);
layerOrder = cell(repeatN * size(result,1), 1);
prePost = cell(repeatN * size(result,1), 1);

for seq_idx = 1:size(result,1)
    seqname = result(seq_idx,1).sequence_name{:};
    
    for layer_idx = 1:size(layers,1)
        if contains(seqname, cell2mat(layers(layer_idx)))
            layerOrder(repeatN * seq_idx - 1) = layers(layer_idx);
            layerOrder(repeatN * seq_idx) = layers(layer_idx);
            break;
        end
        if layer_idx == size(layers,1)
            disp("ERROR: TYPE DOES NOT EXISTS"); break;
        end
    end
    
    if size(seqname) < 10
        len_dirname = size(seqname);
    else
        len_dirname = 10;
    end
    for type_idx = 1:size(types,1)
        if contains(seqname(1:len_dirname), cell2mat(types(type_idx)))
            protein(repeatN * seq_idx - 1) = types(type_idx);
            protein(repeatN * seq_idx) = types(type_idx);
            break;
        end
        if type_idx == size(types,1)
            disp("ERROR: TYPE DOES NOT EXISTS"); break;
        end
    end
    
    prePost(repeatN * seq_idx - 1) = {'pre'};
    prePost(repeatN * seq_idx) = {'post'};

    % volume, pre
    if iscell(result{seq_idx, 2})
        volume(repeatN * seq_idx - 1) = str2double(result{seq_idx, 2}{:}); % pre
    else
        volume(repeatN * seq_idx - 1) = result{seq_idx, 2}; % pre
    end
        
    % volume, post
    if iscell(result{seq_idx, 3})
        volume(repeatN * seq_idx) = str2double(result{seq_idx, 3}{:}); % post
    else
        volume(repeatN * seq_idx) = result{seq_idx, 3}; % post
    end
    
    % SNR, pre
    if iscell(result{seq_idx, 8})
        SNR(repeatN * seq_idx - 1) = str2double(result{seq_idx, 8}{:}); % pre
    else
        SNR(repeatN * seq_idx - 1) = result{seq_idx, 8}; % pre
    end
    
    % SNR, post 
    if iscell(result{seq_idx, 9})
        SNR(repeatN * seq_idx) = str2double(result{seq_idx, 9}{:}); % post
    else
        SNR(repeatN * seq_idx) = result{seq_idx, 9}; % post
    end
end

end