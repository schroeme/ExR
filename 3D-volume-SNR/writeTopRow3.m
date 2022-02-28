function writeTopRow(filename, sheetname)

writematrix('sequence_name', filename, 'Sheet', sheetname, 'Range', 'B1');
writematrix('volume_ref', filename, 'Sheet', sheetname, 'Range', 'C1');
writematrix('volume_pre', filename, 'Sheet', sheetname, 'Range', 'D1');
writematrix('volume_post', filename, 'Sheet', sheetname, 'Range', 'E1');
writematrix('SNR_power_pre', filename, 'Sheet', sheetname, 'Range', 'F1');
writematrix('SNR_power_post', filename, 'Sheet', sheetname, 'Range', 'G1');
writematrix('SNR_pre', filename, 'Sheet', sheetname, 'Range', 'H1');
writematrix('SNR_post', filename, 'Sheet', sheetname, 'Range', 'I1');
writematrix('SNR_power_db_pre', filename, 'Sheet', sheetname, 'Range', 'J1');
writematrix('SNR_power_db_post', filename, 'Sheet', sheetname, 'Range', 'K1');
writematrix('avg_signal_pre', filename, 'Sheet', sheetname, 'Range', 'L1');
writematrix('avg_signal_post', filename, 'Sheet', sheetname, 'Range', 'M1');
writematrix('intersec_pre', filename, 'Sheet', sheetname, 'Range', 'N1');
writematrix('intersec_post', filename, 'Sheet', sheetname, 'Range', 'O1');

end
