%% ASL Study
%
% Inputs:
%   NOISEID  - Path to the .cfl file for the noise scan, generated using
%              BART (v0.0.9) via the TWIXREAD command.
%   DATAID   - Path to the .cfl file for the ASL imaging scan (label and control),
%              generated using BART (v0.0.9) via the TWIXREAD command.
%
% Requirements:
%   - BART Toolbox (v0.0.9)
%     https://github.com/mrirecon/bart
%   - GRAPPA reconstruction tools:
%     https://github.com/mchiew/grappa-tools

noiseDataRaw = readcfl(noiseID);
noiseData = squeeze(noiseDataRaw); % [320 42 160 1 29]
NRO = 2476; NSLICE=71;
[~, NSPIRAL, NCOIL, NVOL] = size(noiseData);
noiseData = reshape(noiseData,NRO,NSLICE,NSPIRAL,NVOL,NCOIL);
EMPTY_IND = find(sum(abs(noiseData(:,:,1,1,1)),1)==0);
SLICE_EMPTY = length(EMPTY_IND);
noiseData = noiseData(:,setdiff(1:NSLICE,EMPTY_IND),:,:,:);
noiseData = reshape(noiseData, [], NSPIRAL, NCOIL, NVOL);
noiseData = permute(noiseData, [1 2 4 3]);
noiseDataReshape = reshape(noiseData, [NRO*(NSLICE-SLICE_EMPTY)*NSPIRAL*NVOL NCOIL]);
T = prewhiten(noiseDataReshape);

dataRaw = readcfl(dataID);

data = squeeze(dataRaw);
NRO = 2476; NSLICE=71;
[~, NSPIRAL, NCOIL, NVOL] = size(data);
data = permute(data, [1 2 4 3]);
data = reshape(data, NRO, NSLICE, NSPIRAL, NVOL, NCOIL); 
d = data(:,EMPTY_IND,:,:,:);
data = data(:, setdiff(1:NSLICE,EMPTY_IND), :, :, :);


NVOL = size(data,4);

ACTUAL_NSLICE = NSLICE - SLICE_EMPTY;
getSliceIdx = @(i) lookupTable(1:ACTUAL_NSLICE,1:ACTUAL_NSLICE,i);

KEY_TBL = [];
xStep=40;zStep=10;
for k=1:zStep:ACTUAL_NSLICE
    zRange = k:min([k+zStep-1,ACTUAL_NSLICE]);
    zIdx = [];
    for i = zRange
        zIdx(end+1) = getSliceIdx(i);
    end
    for j = 1:NSPIRAL
        for i = 1:xStep:NRO
        xRange = i:min([i+xStep-1,NRO]);
        KEY_TBL = unique([KEY_TBL length(xRange)*length(zRange)*NCOIL]);
        end
    end
end

N = 1000;
SMC_TBL=[];
for x = KEY_TBL
    y=NVOL;
    SMC = zeros(N,1);
    for i = 1:N
    [~,tmp,~] = svd(...
     randn([x y])*sqrt(0.5) + randn([x y])*sqrt(0.5)*1j,'econ');
    SMC(i) = tmp(1,1);
    end
    SMC_TBL = [SMC_TBL mean(SMC)];
end


getThreshold = @(i) lookupTable(KEY_TBL,SMC_TBL,i);


dataDe = zeros(size(data));
pct_remove = zeros(length(1:xStep:NRO),length(1:NSPIRAL),length(1:zStep:ACTUAL_NSLICE));
for k=1:zStep:ACTUAL_NSLICE
    zRange = k:min([k+zStep-1,ACTUAL_NSLICE]);
    zIdx = [];
    for i = zRange
        zIdx(end+1) = getSliceIdx(i);
    end
    dataPrewhiten = reshape(data(:,zIdx,:,:,:), [], NCOIL)*T;
    dataPrewhiten = reshape(dataPrewhiten, [NRO, length(zRange), NSPIRAL, NVOL, NCOIL]);
    removeCounts = 0;
    fullCounts = 0;
    for j = 1:NSPIRAL
        for i = 1:xStep:NRO
            xRange = i:min([i+xStep-1,NRO]);
            casorati = dataPrewhiten(xRange,1:length(zRange),j,1:NVOL,1:NCOIL);
            casorati = reshape(permute(casorati, [1 2 3 5 4]),[],NVOL);
            [U,S,V] = svd(casorati,'econ');
            SMC = getThreshold(size(casorati,1));
            removeCounts = removeCounts+sum(diag(S)<SMC);
            fullCounts = fullCounts + size(S,1);
            S(S<SMC) = 0;
            casoratiDe = U*S*V';
            casoratiDe = reshape(casoratiDe,length(xRange),length(zRange),1,NCOIL, NVOL);
            casoratiDe = permute(casoratiDe,[1 2 3 5 4]);
            pct_remove(ceil(i/xStep),j,ceil(k/zStep)) = sum(diag(S)<SMC)/size(S,1);
            casoratiDe = reshape(casoratiDe,[],NCOIL)*inv(T);
            casoratiDe = reshape(casoratiDe,[length(xRange), length(zRange), 1, NVOL, NCOIL]);
            dataDe(xRange, zRange, j, 1:NVOL, 1:NCOIL) = casoratiDe;
        end
    end
    fprintf('%d out of %d finished - %.2f %% removed.\n', zRange(end),ACTUAL_NSLICE,removeCounts/fullCounts*100);
end
fulldata = zeros(NRO, NSLICE, NSPIRAL, NVOL, NCOIL); 
fulldata(:,setdiff(1:NSLICE,EMPTY_IND),:,:,:) = dataDe; % Denoised k-space data