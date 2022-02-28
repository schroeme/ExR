function im_norm = getNormImg(im)

im_norm = im - min(im, [], 'all');
im_norm = im_norm / max(im_norm, [], 'all');

end