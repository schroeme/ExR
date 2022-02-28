function scores = getScore(pre, pre_result, post, post_result, ref_result, STD_pre, STD_post)

% compute volume
volume_ref = sum(ref_result, 'all'); % score at 1
volume_pre = sum(pre_result, 'all'); % score at 2
volume_post = sum(post_result, 'all'); % score at 3

% compute signal and noise
signal_pre = pre(pre_result);
signal_post = post(post_result);
% noise_pre = pre(~ pre_result);
% noise_post = post(~ post_result);

% compute signal power and noise power
% signal_power_pre = signal_pre .^ 2;
% signal_power_post = signal_post .^ 2;
% noise_power_pre = noise_pre .^ 2;
% noise_power_post = noise_post .^ 2;

% compute SNR based on power and NOT power
% SNR_power_pre = mean(signal_power_pre, 'all') / mean(noise_power_pre, 'all'); % score at 4
% SNR_power_post = mean(signal_power_post, 'all') / mean(noise_power_post, 'all'); % score at 5
SNR_pre = mean(signal_pre, 'all') / STD_pre; % score at 6
SNR_post = mean(signal_post, 'all') / STD_post; % score at 7
% SNR_power_db_pre = 10 * log10 (SNR_power_pre); % score at 8
% SNR_power_db_post = 10 * log10 (SNR_power_post); % score at 9

% compute intensity
% avg_signal_pre = mean(signal_pre, 'all'); % score at 10
% avg_signal_post = mean(signal_post, 'all'); % score at 11

% compute intersection
% intersec_pre = sum(pre_result .* ref_result, 'all') / volume_ref; % score at 12
% intersec_post = sum(post_result .* ref_result, 'all') / volume_ref; % score at 13

scores = [volume_ref, volume_pre, volume_post, ...
    SNR_pre, SNR_post];

end