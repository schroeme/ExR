function [STD_pre, STD_post] = getSTD(STD, type_order, layer_order, ...
    AB_order, seqname, layername)

if size(seqname) < 10
    len_dirname = size(seqname);
else
    len_dirname = 10;
end
for type_idx = 1:size(type_order,1)
    if contains(seqname(1:len_dirname), cell2mat(type_order(type_idx)))
        column_idx = 2 * (type_idx - 1) + 1;
        break;
    end
    if type_idx == size(type_order,1)
        disp("ERROR: TYPE DOES NOT EXISTS"); break;
    end
end

for layer_idx = 1:size(layer_order,1)
    if contains(layername, cell2mat(layer_order(layer_idx)))
        row_idx = layer_idx;
        break;
    end
    if layer_idx == size(layer_order,1)
        disp("ERROR: TYPE DOES NOT EXISTS"); break;
    end
end

for AB_idx = 1:size(AB_order,1)
    if contains(seqname(1), cell2mat(AB_order(AB_idx)))
        row_idx = row_idx + size(layer_order,1) * (AB_idx - 1);
        break;
    end
    if AB_idx == size(AB_order,1)
        disp("ERROR: TYPE DOES NOT EXISTS"); break;
    end
end

STD_pre = STD(row_idx, column_idx);
STD_post = STD(row_idx, column_idx + 1);

end