function [n_channel, n_z_stack, n_h, n_w] = getImgInfo(imname)

% get number of channels and stacks
iminfo = imfinfo(imname);
iminfo_detail = iminfo(1).ImageDescription;
iminfo_detail = splitlines(iminfo_detail);
for info_idx = 1:size(iminfo_detail, 1)
    iminfo_split = split(iminfo_detail{info_idx, 1}, '=');
    if strcmp(iminfo_split{1}, 'channels')
        n_channel = str2double(iminfo_split{2});
    end
    if strcmp(iminfo_split{1}, 'slices')
        n_z_stack = str2double(iminfo_split{2});
    end
end

n_h = iminfo(1).Height;
n_w = iminfo(1).Width;

if size(iminfo, 1) ~= n_channel * n_z_stack
    disp ("ERROR: number of channels does not match"); return;
end

end