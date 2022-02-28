close all; clear all; clc;

outputFile = 'result.xlsx';
resultFile = dir('result/*.xlsx');
sheetname_set = {'n=7'};

% iterate over sheets
for sheet_idx = 1:size(sheetname_set,1)
    sheetname = sheetname_set{sheet_idx};
    
    % iterate over result files
    for file_idx = 1:size(resultFile,1)
        filename = resultFile(file_idx).name;
       
        writeTopRow4(outputFile, sheetname);
        
        result = readtable(fullfile('result', filename), 'Sheet', sheetname);
        avg_result = result(size(result,1)-4,:);
        std_result = result(size(result,1)-2,:);
        se_result = result(size(result,1),:);
        
        writematrix(filename(1:end-12), outputFile, 'Sheet', sheetname, ...
            'Range', sprintf('A%d', file_idx+1));
        writetable(avg_result(1, 2:3), outputFile, 'Sheet', sheetname, ...
            'Range', sprintf('B%d', file_idx+1), 'WriteVariableNames', 0);
        writetable(std_result(1, 2:3), outputFile, 'Sheet', sheetname, ...
            'Range', sprintf('D%d', file_idx+1), 'WriteVariableNames', 0);
        writetable(se_result(1, 2:3), outputFile, 'Sheet', sheetname, ...
            'Range', sprintf('F%d', file_idx+1), 'WriteVariableNames', 0);
        writetable(avg_result(1, 4:5), outputFile, 'Sheet', sheetname, ...
            'Range', sprintf('H%d', file_idx+1), 'WriteVariableNames', 0);
        writetable(std_result(1, 4:5), outputFile, 'Sheet', sheetname, ...
            'Range', sprintf('J%d', file_idx+1), 'WriteVariableNames', 0);
        writetable(se_result(1, 4:5), outputFile, 'Sheet', sheetname, ...
            'Range', sprintf('L%d', file_idx+1), 'WriteVariableNames', 0);
    end
end