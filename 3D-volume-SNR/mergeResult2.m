close all; clear all; clc;

types = {'Ca cha'; 'PSD95'; 'syngap'; 'homer'; 'Rim'; 'basoon'; 'shank'};
layers = {'L1'; 'L23'; 'L4'};
AB = {'A'; 'B'};

outputFile = 'result_final.xlsx';
filename = 'result.xlsx';
sheetname_set = {'n=7'};

% iterate over sheets
for sheet_idx = 1:size(sheetname_set,1)
    sheetname = sheetname_set{sheet_idx};
    
    writeTopRow5(outputFile, sheetname);
    
    result = readtable(filename, 'Sheet', sheetname);
    
    avg_result = zeros(size(AB,1) * size(types,1), 8);
    for seq_idx = 1:size(result,1)
        seqname = result(seq_idx,1).sequence_name{:};
        
        for AB_idx = 1:size(AB,1)
            if contains(seqname(1), cell2mat(AB(AB_idx)))
                row_idx = size(types,1) * (AB_idx - 1); break;
            end
            if AB_idx == size(AB,1)
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
                row_idx = row_idx + type_idx; break;
            end
            if type_idx == size(types,1)
                disp("ERROR: TYPE DOES NOT EXISTS"); break;
            end
        end

        for score_idx = 2:size(result,2)
            if iscell(result{seq_idx, score_idx})
                avg_result(row_idx, score_idx-1) = avg_result(row_idx, score_idx-1) + ...
                    str2double(result{seq_idx, score_idx}{:}) / size(layers,1);
            else
                avg_result(row_idx, score_idx-1) = avg_result(row_idx, score_idx-1) + ...
                    result{seq_idx, score_idx} / size(layers,1);
            end
        end
    end
    
    writematrix(repmat(AB{1}, size(types,1),1), outputFile, 'Sheet', sheetname, ...
        'Range', sprintf('A%d', 2));
    writematrix(repmat(AB{2}, size(types,1),1), outputFile, 'Sheet', sheetname, ...
        'Range', sprintf('A%d', size(types,1)+2));
    writecell(repmat(types,2,1), outputFile, 'Sheet', sheetname, ...
        'Range', sprintf('B%d', 2));
    writematrix(avg_result, outputFile, 'Sheet', sheetname, ...
        'Range', sprintf('C%d', 2));
end



