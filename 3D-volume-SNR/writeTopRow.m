function writeTopRow(filename, sheetname)

writematrix('volume_ref', filename, 'Sheet', sheetname, 'Range', 'A1');
writematrix('volume_pre', filename, 'Sheet', sheetname, 'Range', 'B1');
writematrix('volume_post', filename, 'Sheet', sheetname, 'Range', 'C1');
writematrix('SNR_pre', filename, 'Sheet', sheetname, 'Range', 'D1');
writematrix('SNR_post', filename, 'Sheet', sheetname, 'Range', 'E1');

end