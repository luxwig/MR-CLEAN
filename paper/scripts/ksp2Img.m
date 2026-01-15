function imgData = ksp2Img(ksp, dims)
    imgData = ksp;
    for i = dims
        imgData = ifftshift(ifft(imgData,[],i),i)*size(imgData,i);
    end
end