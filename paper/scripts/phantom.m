%% Phantom Study
%
% Inputs:
%   NOISEID  - Path to the TWIX file for the noise scan.
%   DATAID   - Path to the TWIX file for the imaging scan.
%
% Requirements:
%   - mapVBVD
%     https://github.com/pehses/mapVBVD

noiseRaw = mapVBVD(noiseID);

noiseData = squeeze(noiseRaw{2}.image(:,:,:,1,:,1,1,1,:,1,:)); 
noiseData = sum(noiseData,6);
noiseData = reshape(noiseData,size(noiseData,1),size(noiseData,2),size(noiseData,3),[]);
noiseData = permute(noiseData, [3 4 1 2]); % PE x slice x RO x Coil
[NPE, NSLICE, NRO, NCOIL] = size(noiseData);
noiseData = reshape(noiseData, [NPE*NSLICE NRO NCOIL]);
noiseFreq = fftshift(fft(ifftshift(noiseData,2),[],2)/NRO,2);
noiseFreq = noiseFreq(:,NRO/4+1:NRO/4+NRO/2,:);
noiseCorrect = ifftshift(ifft(fftshift(noiseFreq,2),[],2)*NRO/2,2);
noiseCorrectReshape = reshape(noiseCorrect, [NPE*NSLICE*NRO/2 NCOIL]);
T = prewhiten(noiseCorrectReshape);


dataRaw = mapVBVD(dataID);


xStep = 6; yStep = 6; zStep = 6;
NRO = dataRaw{2}.image.NCol;
NPE = dataRaw{2}.image.NLin;
NSLICE = dataRaw{2}.image.NPar;
NVOL = dataRaw{2}.image.NRep;
NVOL = 35;
NCOIL = dataRaw{2}.image.NCha;

SMC_IND = [];
for k=1:zStep:NSLICE
    zRange = k:min([k+zStep-1,NSLICE]);
    for i = 1:xStep:NRO/2
        xRange = i:min([i+xStep-1,NRO/2]);
        for j = 1:yStep:NPE
            yRange = j:min([j+yStep-1,NPE]);
            SMC_IND = unique([SMC_IND length(xRange)*length(yRange)*length(zRange)*NCOIL]);
        end
    end
end


N = 1000;
SMC_TBL=[];
KEY_TBL=SMC_IND;
for x = SMC_IND
    y=NVOL;
    SMC = zeros(N,1);
    for i = 1:N
    [~,tmp,~] = svd(...
    ... % randn([x y])*sqrt(0.5) + randn([x y])*sqrt(0.5)*1j,'econ');
    randn([x y]),'econ');
    SMC(i) = tmp(1,1);
    end
    SMC_TBL = [SMC_TBL mean(SMC)];
end

getSliceIdx = @(i) lookupTable(1:NSLICE,1:NSLICE,i);
getThreshold = @(i) lookupTable(KEY_TBL,SMC_TBL,i);

d = zeros(NRO/2,NPE,NSLICE,NVOL,NCOIL);
for k=1:zStep:NSLICE
    zRange = k:min([k+zStep-1,NSLICE]);
    zIdx = [];
    for i = zRange
        zIdx(end+1) = getSliceIdx(i);
    end
    data = squeeze(dataRaw{2}.image(:,:,:,zIdx,1,1,1,1,1:NVOL));
    data = permute(data, [1 3 4 5 2]);
    dataPrewhiten = fftshift(fft(ifftshift(data,1),[],1)/NRO,1);
    dataPrewhiten = dataPrewhiten(NRO/4+1:NRO/4+NRO/2,:,:,:,:);
    dataPrewhiten = ifftshift(ifft(fftshift(dataPrewhiten,1),[],1)*NRO/2,1);
    dataPrewhiten = reshape(dataPrewhiten, [], NCOIL)*T;
    dataPrewhiten = reshape(dataPrewhiten, [NRO/2, NPE, length(zRange), NVOL, NCOIL]);
    dataDe = zeros(size(dataPrewhiten));
    f = waitbar(0,sprintf('Slice %d - %d, 0%%', zRange(1), zRange(end)));
    pg = 0;
    removeCounts = 0;
    fullCounts = 0;
    for i = 1:xStep:NRO/2
        xRange = i:min([i+xStep-1,NRO/2]);
        for j = 1:yStep:NPE
            yRange = j:min([j+yStep-1,NPE]);
            pg = pg + length(yRange)*length(xRange);
            casorati = dataPrewhiten(xRange,yRange,1:length(zRange),1:NVOL,1:NCOIL);
            casorati = reshape(permute(casorati, [1 2 3 5 4]),[],NVOL);
            [U,S,V] = svd(casorati, 'econ');
            SMC = getThreshold(size(casorati,1));
            removeCounts = removeCounts+sum(diag(S)<SMC);
            fullCounts = fullCounts + size(S,1);
            S(S<SMC) = 0;
            casoratiDe = U*S*V';
            casoratiDe = reshape(casoratiDe,length(xRange),length(yRange),length(zRange),NCOIL, NVOL);
            casoratiDe = permute(casoratiDe,[1 2 3 5 4]);
            dataDe(xRange, yRange, 1:length(zRange), 1:NVOL, 1:NCOIL) = casoratiDe;
            f = waitbar(pg/(NRO/2*NPE),f,sprintf('Slice %d - %d, %.2f%%', zRange(1), zRange(end),pg/(NRO/2*NPE)*100));
        end
    end
    close(f);
    d(:,:,zRange,:,:) = dataDe;
    fprintf('%d out of %d finished - %.2f %% removed.\n', zRange(end),NSLICE,removeCounts/fullCounts*100);
end

dataFull_recon = getImage(reshape(reshape(d,[],NCOIL)*inv(T),NRO/2,NPE,NSLICE,NVOL,NCOIL)); % Denoised image in image space