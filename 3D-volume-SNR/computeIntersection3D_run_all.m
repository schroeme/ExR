close all; clear all; clc;

% order in STD file
type_order = {'cacha'; 'psd95'; 'syngap'; 'homer'; 'rim'; 'bassoon'; 'shank'};
layer_order = {'L1'; 'L23'; 'L4'};
AB_order = {'A'; 'B';'C'};

% read STD from excel file % IMPORTANT --- keep table orders ---
STDfile = readmatrix('./STD_revision.xlsx');
STD = STDfile(:,3:end);

%% read sequence
for thresh_scale = 7:7
    seqname_set = dir('/Users/margaretschroeder/Dropbox (MIT)/BoydenLab/idExMDecrowdingSegs/Decrowding crop image/images');
    parfor seq_idx = 4:size(seqname_set, 1) % temporary to exclude couple of dir
        seqname = seqname_set(seq_idx, 1).name;
        disp(seqname);
        
        if contains(seqname, 'A1_cacha_shank') | ...
                contains(seqname, 'B1_cacha_shank') | ...
                contains(seqname, 'A3_syngap_shank') | ...
                contains(seqname, 'B3_syngap_shank') | ...
                contains(seqname, 'A4_homer_shank') | ...
                contains(seqname, 'B4_homer_shank')
                
            pre_idx = 2; post_idx = 3;
        else
            pre_idx = 3; post_idx = 2;
%         else
%             disp("ERROR: check file name");
        end
        
        dirname = fullfile('/Users/margaretschroeder/Dropbox (MIT)/BoydenLab/idExMDecrowdingSegs/Decrowding crop image/images', seqname);
        layername_set = dir(fullfile(dirname, '*L*'));
        
        %% iterate layers
        for layer_idx = 1:size(layername_set, 1)
            layername = layername_set(layer_idx, 1).name;
            curr_dir_name = fullfile(dirname, layername);
            disp (layername); 
            
            % write top row to excel file
            filename = sprintf('%s_%s_result.xlsx', seqname, layername);
            sheetname = sprintf('n=%d', thresh_scale);
            writeTopRow(filename, sheetname);
            
            % read tif file names
            files = dir(fullfile(curr_dir_name, '*.tif'));
            if size(files, 1) == 0
                disp ("ERROR: no tif file exists in "); disp(curr_dir_name);
            end
            
            % get STD values
            [STD_pre, STD_post] = getSTD(STD, type_order, layer_order, AB_order, seqname, layername);
            
            %% iterate images
%             avg_score = zeros(1,5);
            cumu_score = zeros(size(files, 1), 5);
            for im_idx = 1:size(files, 1)
                imname = fullfile(curr_dir_name, files(im_idx).name);
                
                % get number of channels and stacks
                [n_channel, n_z_stack, n_h, n_w] = getImgInfo(imname);
                
                % read image % channel then, stack % IMPORTANT --- image order ---
                ref = zeros(n_h, n_w, n_z_stack);
                pre = zeros(n_h, n_w, n_z_stack);
                post = zeros(n_h, n_w, n_z_stack);
                for idx_stack = 1: n_z_stack
                    ref(:,:,idx_stack) = imread(imname, 'tiff', n_channel * (idx_stack - 1) + 1);
                    pre(:,:,idx_stack) = imread(imname, 'tiff', n_channel * (idx_stack - 1) + pre_idx);
                    post(:,:,idx_stack) = imread(imname, 'tiff', n_channel * (idx_stack - 1) + post_idx);
                end
                
                % subtract background signal
                ref_norm = ref - min(ref, [], 'all');
                pre_norm = pre - min(pre, [], 'all');
                post_norm = post - min(post, [], 'all');
                
                % determine threshold % ------- IMPORTANT ------- *
                T_ref = thresh_scale * ((STD_pre + STD_post) / 2) ;
                T_pre = thresh_scale * STD_pre;
                T_post = thresh_scale * STD_post;
                
                ref_result = connectedComponent3D(ref_norm, T_ref);
                pre_result = connectedComponent3D(pre_norm, T_pre);
                post_result = connectedComponent3D(post_norm, T_post);
                
%                 % write images
%                 curr_write_dir = fullfile('result_fig', seqname, layername, num2str(thresh_scale));
%                 if ~isfolder(curr_write_dir)
%                     mkdir(curr_write_dir)
%                 end
%                 for z_idx = 1:n_z_stack
%                     scale = max([ref; pre; post], [], 'all');
%                     imwrite(ind2rgb(ref(:,:,z_idx), hot(scale)), fullfile(curr_write_dir, sprintf('%s_z%d_ref_orig.png', files(im_idx).name(1:end-4), z_idx)));
%                     imwrite(ind2rgb(pre(:,:,z_idx), hot(scale)), fullfile(curr_write_dir, sprintf('%s_z%d_pre_orig.png', files(im_idx).name(1:end-4), z_idx)));
%                     imwrite(ind2rgb(post(:,:,z_idx), hot(scale)), fullfile(curr_write_dir, sprintf('%s_z%d_post_orig.png', files(im_idx).name(1:end-4), z_idx)));
%                     imwrite(ref_result(:,:,z_idx), fullfile(curr_write_dir, sprintf('%s_z%d_ref_result.png', files(im_idx).name(1:end-4), z_idx)));
%                     imwrite(pre_result(:,:,z_idx), fullfile(curr_write_dir, sprintf('%s_z%d_pre_result.png', files(im_idx).name(1:end-4), z_idx)));
%                     imwrite(post_result(:,:,z_idx), fullfile(curr_write_dir, sprintf('%s_z%d_post_result.png', files(im_idx).name(1:end-4), z_idx)));
%                 end
                
                % compute score
                scores = getScore(pre, pre_result, post, post_result, ref_result, STD_pre, STD_post);
                
                % write result to excel
                writematrix(scores, filename, 'Sheet', sheetname, ...
                    'Range', sprintf('A%d', im_idx+1));
                
                % compute average score
                cumu_score(im_idx, :) = scores;
%                 avg_score = avg_score + scores / size(files,1);
                
                if im_idx == size(files,1)
                    avg_score = mean(cumu_score);
                    writematrix('mean', filename, 'Sheet', sheetname, 'Range', sprintf('A%d', im_idx+2));
                    writematrix(avg_score, filename, 'Sheet', sheetname, ...
                        'Range', sprintf('A%d', im_idx+3));
                    
                    std_score = std(cumu_score);
                    writematrix('std', filename, 'Sheet', sheetname, 'Range', sprintf('A%d', im_idx+4));
                    writematrix(std_score, filename, 'Sheet', sheetname, ...
                        'Range', sprintf('A%d', im_idx+5));
                    
                    se_score = std_score / sqrt(size(cumu_score,1));
                    writematrix('se', filename, 'Sheet', sheetname, 'Range', sprintf('A%d', im_idx+6));
                    writematrix(se_score, filename, 'Sheet', sheetname, ...
                        'Range', sprintf('A%d', im_idx+7));
                    
                end
                
                % ============ visualize ============
                %         visualizeResult(ref, pre_result, post_result);
                % ============ visualize ============
            end
        end
    end
end
disp('DONE');
