function writeTopRow(filename, sheetname)

writematrix('sequence_name', filename, 'Sheet', sheetname, 'Range', 'A1');
writematrix('volume_pre', filename, 'Sheet', sheetname, 'Range', 'B1');
writematrix('volume_post', filename, 'Sheet', sheetname, 'Range', 'C1');
writematrix('volume_pre_STD', filename, 'Sheet', sheetname, 'Range', 'D1');
writematrix('volume_post_STD', filename, 'Sheet', sheetname, 'Range', 'E1');
writematrix('volume_pre_SE', filename, 'Sheet', sheetname, 'Range', 'F1');
writematrix('volume_post_SE', filename, 'Sheet', sheetname, 'Range', 'G1');
writematrix('SNR_pre', filename, 'Sheet', sheetname, 'Range', 'H1');
writematrix('SNR_post', filename, 'Sheet', sheetname, 'Range', 'I1');
writematrix('SNR_pre_STD', filename, 'Sheet', sheetname, 'Range', 'J1');
writematrix('SNR_post_STD', filename, 'Sheet', sheetname, 'Range', 'K1');
writematrix('SNR_pre_SE', filename, 'Sheet', sheetname, 'Range', 'L1');
writematrix('SNR_post_SE', filename, 'Sheet', sheetname, 'Range', 'M1');

end
