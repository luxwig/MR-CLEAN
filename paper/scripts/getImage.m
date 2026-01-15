function dataTobeWrite = getImage(kspaceData)
    dataTobeWrite = kspaceData;
    dataDim = [1 2 3];

    dataTobeWrite = squeeze(dataTobeWrite);
    for i = dataDim
        dataTobeWrite = ifftshift(ifft(fftshift(dataTobeWrite,i),[],i),i);
    end
    dataTobeWrite = sqrt(sum(abs(dataTobeWrite).^2,5));
end