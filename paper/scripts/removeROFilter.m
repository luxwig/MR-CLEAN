function data = removeROFilter(data)
    NRO = size(data,1);
    data = fftshift(fft(ifftshift(data,1),[],1)/NRO,1);
    data = data(NRO/4+1:NRO/4+NRO/2,:,:,:,:);
    data = ifftshift(ifft(fftshift(data,1),[],1)*NRO/2,1);
end
