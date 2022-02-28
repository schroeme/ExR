function visualizeResult(ref, pre_result, post_result)

figure(1); volshow(ref); movegui('north');
% volumeViewer(ref);
figure(2); sliceViewer(ref, 'ScaleFactors', [10,10,1]); movegui('center');
figure(3); labelvolshow(pre_result); movegui('northwest');
figure(4); sliceViewer(pre_result, 'ScaleFactors', [10,10,1]); movegui('west');
figure(5); labelvolshow(post_result); movegui('northeast');
figure(6); sliceViewer(post_result, 'ScaleFactors', [10,10,1]); movegui('east');

end