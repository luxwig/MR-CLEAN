%% Cardiac Study
%
% Inputs:
%   NOISEID  - Path to the .cfl file for the noise scan,
%              generated using BART (v0.0.9) via the TWIXREAD command.
%   DATAID   - Path to the .cfl file for the imaging scan,
%              generated using BART (v0.0.9) via the TWIXREAD command.
%
% Requirements:
%   - BART Toolbox (v0.0.9)
%     https://github.com/mrirecon/bart
%   - GRAPPA reconstruction tools:
%     https://github.com/mchiew/grappa-tools

noiseDataRaw = readcfl(noiseID);

noiseData = squeeze(noiseDataRaw);
noiseData = noiseData(:,:,:,1:23,:);
m = abs(sum(squeeze(noiseData(:,:,1,1,1))))~=0;
noiseData = noiseData(:,m,:,:,:);
noiseData = permute(noiseData,[1 3 2 4 5]);
[NRO NCOIL] = size(noiseData,[1 2]);
noiseData = reshape(noiseData, NRO, NCOIL, []);
noiseCorrect = removeROFilter(noiseData);
noiseCorrect = permute(noiseCorrect, [3 1 2]);
noiseCorrectReshape = reshape(noiseCorrect, [], NCOIL);
T = prewhiten(noiseCorrectReshape);

dataRaw = readcfl(dataID);

NRO_FULL = 240;
NPE_FULL = 176;
CALIB_LINES = [65:108]+1;
NVOL = 23;
data = squeeze(dataRaw);
data = data(:,:,:,1:NVOL,:);
IND = abs(sum(squeeze(data(:,:,1,1,1))))~=0;
data = data(:,IND,:,:,:);
IND = [0 IND];
[NRO NPE NCOIL NSLICE] = size(data, [1 2 3 5]);
data = permute(data, [1 2 5 4 3]); % RO PE SLICE VOL COIL

xStep = 4; yStep = 4; zStep = 1;

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
    [~,tmp,~] = svd(randn([x y])*sqrt(0.5) + randn([x y])*sqrt(0.5)*1j,'econ');
    SMC(i) = tmp(1,1);
    end
    SMC_TBL = [SMC_TBL mean(SMC)];
end

getSliceIdx = @(i) lookupTable(1:NSLICE,1:NSLICE,i);
getThreshold = @(i) lookupTable(KEY_TBL,SMC_TBL,i);

pct_remove = zeros(length(1:xStep:NRO/2),length(1:yStep:NPE),length(1:zStep:NSLICE));
dataPrewhiten = removeROFilter(data);
dataPrewhiten = reshape(dataPrewhiten, [], NCOIL)*T;
dataPrewhiten = reshape(dataPrewhiten, [NRO/2, NPE, NSLICE, NVOL, NCOIL]);
dataDe = zeros(size(dataPrewhiten));

for k=1:zStep:NSLICE
    zRange = k:min([k+zStep-1,NSLICE]);
    zIdx = [];
    for i = zRange
        zIdx(end+1) = getSliceIdx(i);
    end
    
    removeCounts = 0;
    fullCounts = 0;
    for i = 1:xStep:NRO/2
        xRange = i:min([i+xStep-1,NRO/2]);
        for j = 1:yStep:NPE
            yRange = j:min([j+yStep-1,NPE]);
            casorati = dataPrewhiten(xRange,yRange,zRange,1:NVOL,1:NCOIL);
            casorati = reshape(permute(casorati, [1 2 3 5 4]),[],NVOL);
            [U,S,V] = svd(casorati,'econ');
            SMC = getThreshold(size(casorati,1));
            removeCounts = removeCounts+sum(diag(S)<SMC);
            fullCounts = fullCounts + size(S,1);
            S(S<SMC) = 0;
            casoratiDe = U*S*V';
            casoratiDe = reshape(casoratiDe,length(xRange),length(yRange),length(zRange),NCOIL, NVOL);
            casoratiDe = permute(casoratiDe,[1 2 3 5 4]);
            pct_remove(ceil(i/xStep),ceil(j/yStep),ceil(k/zStep)) = sum(diag(S)<SMC)/size(S,1);
            dataDe(xRange, yRange, zRange, 1:NVOL, 1:NCOIL) = casoratiDe;
        end
    end
    fprintf('%d out of %d finished - %.2f %% removed.\n', zRange(end),NSLICE,removeCounts/fullCounts*100);
end

dataDe = reshape(dataDe, [], NCOIL)*inv(T);
dataDe = reshape(dataDe, [NRO/2, NPE, NSLICE, NVOL, NCOIL]);
dataFull = zeros(NRO/2,NPE_FULL,NSLICE, NVOL, NCOIL);
dataFull(:,IND==1,:,:,:) = dataDe;
dataFull = permute(dataFull,[5 1 2 3 4]);

dataFull_grappa = zeros(size(dataFull));
tmp = zeros(size(dataFull_grappa(:,:,:,1,:)));
for i = 1:size(dataFull_grappa,4)
    tmp(:,:,(1:4:NPE_FULL)+3,:,:) = dataFull(:,:,(1:4:NPE_FULL)+3,i,:);
    dataFull_grappa(:,:,:,i,:)= ...
                 grappa(tmp, ...
                 dataFull(:,:,CALIB_LINES,i,:),[1 4],[4 4]);
    fprintf('GRAPPA - %d out of %d finished.\n',i,size(dataFull_grappa,4));
end
dataFull_grappa(:,:,CALIB_LINES,:,:) = dataFull(:,:,CALIB_LINES,:,:);

dataFull_recon = size(dataFull_grappa);
dataFull_recon(2) = NRO_FULL;
dataFull_recon = zeros(dataFull_recon);
dataFull_recon(:,end-size(dataFull_grappa,2)+1:end,:,:,:) = dataFull_grappa;
dataFull_recon = ksp2Img(dataFull_recon,[2 3]);
dataFull_recon = sum(abs(dataFull_recon).^2,1);
dataFull_recon = squeeze(dataFull_recon); % Denoised image in image space