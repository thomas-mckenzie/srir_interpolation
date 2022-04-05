% close all
clear
clc

% sweep
irPath = '/Users/mckenzt1/Documents/RT Dataset/rt dataset SOFA/SOFA Files/Meeting Room to Hallway SOFA - Dataset of Room Impulse Responses for the Transition between Coupled Rooms/';
irName = 'Room Transition RIRs_Meeting Room to Hallway_Source in Room_No Line of Sight.sofa';

% positions are in the motive system with the "z up" setting
% which is the same, right handed, x in the front system that
% sofa and the tv convolver use
% x front y left z up

sofa = SOFAload([irPath,irName]);

%% configure

% The amount to which to downsample the original measurements
% (2 would mean using 1 of every 2, 10 means 1 in every 10,
% 25 means 1 in every 25, so for 101 measurements it would just be using 5)
measReduction = 10;

% minPhase, rotationOnly or fixedSpectrum or meanSpectrum
INTERPOLATION_MODE = 'meanSpectrum';

directSoundLength_samp = 200; % nils had at 200
fade_samp_DS2ER = 20; % DS direct sound ER early refs LR late reverb
fade_samp_ER2LR = sofa.Data.SamplingRate / 100;
pre_samp = 5; % for truncation to direct sound
pre_samp2 = 0; % 2nd truncation - in NODELAY

earlyRefsCutoff_seconds = 0.2;

% mods, should not be needed
directMod_db = 0;
roomMod_db = 0;

% SH order for interpolated signals
NshInterp = 4;

% FOR CASE OF ROOM TRANSITION (hard coded)
positionsCorrect = round(sofa.ListenerPosition*100);
positionsCorrect(:,3) = 155;
positions = positionsCorrect(1:measReduction:end,:);

%% init

% work out whether interpolation is 1D, 2D or 3D
isPlaneUsed=zeros(length(positions(1,:)),1);
for i = 1:length(positions(1,:))
    if all(positions(:,i) == positions(1,i))
        isPlaneUsed(i) = 0;
    else
        isPlaneUsed(i) = 1;
    end
end
mode_dimension = sum(isPlaneUsed); % 1 = 1D (planar), 2 = 2D, 3 = 3D

% calculate the positions
positionsInterpolatedX_cm = positionsCorrect(:,1);
positionsX_cm = positionsInterpolatedX_cm(1:measReduction:end);
receiverHeight_cm = positionsCorrect(1,3);

if mode_dimension == 1 % hard coded - lazy - needs fix
    positionsInterpolatedY_cm = positionsCorrect(1,2);
    positionsY_cm = positionsCorrect(1,2);
else
    positionsInterpolatedY_cm = positionsCorrect(:,2);
    positionsY_cm = positionsInterpolatedY_cm(1:measReduction:end);
end

hCorrect = permute(sofa.Data.IR,[3,2,1]);
fs = sofa.Data.SamplingRate;
irLength_samp = length(hCorrect(:,1,1));

numMeasurements = length(positionsX_cm) * length(positionsY_cm);
numAmbiChannels = length(hCorrect(1,:,1));
hAll = hCorrect(:,:,1:measReduction:101);

nfft = 2^nextpow2(directSoundLength_samp);

% target=1/exp(1);
target=0.1;

for i = 1:length(hAll(1,1,:))
    [decay,normalisation_value] = rir2decay(hAll(:,1,i), fs, 1000, 1, 1, 1);
    t = (0:length(decay)-1)./fs;
    
    ind_at_t_L(i) = find(decay<target,1);
    % t_L_samp(i) = 2.^nextpow2(ind_at_t_L); % next power of 2
    earlyRefsCutoff_samp(i) = ceil(ind_at_t_L(i)/1000)*1000; % nearest 1000
    
end

% earlyRefsCutoff_samp = earlyRefsCutoff_seconds * sofa.Data.SamplingRate;
% earlyRefsCutoff_samp = earlyRefsCutoff_seconds * sofa.Data.SamplingRate;



%% initial DOA

cutStart_samp = 1; % use this to bypass any additional time alignment



m = SRIR_Measurement();
m.setShMicArray(4, 0.032);
m.fs = fs;

directSoundDirections = zeros(numMeasurements, 3);
aziMeasuredDeg = zeros(numMeasurements, 1);
for iMeas = 1:numMeasurements
    m.srir = hAll(1:directSoundLength_samp, :, iMeas);
    opts.analysisMethod = 'iv';
    a = SRIR_Analysis(m, opts);
    a.run();
    
    idxNonZero = find(sum(a.doa~=0, 2));
    doa = mean(a.doa(idxNonZero, :));
    doa = doa ./ vecnorm(doa, 2, 2);
    
    directSoundDirections(iMeas, :) = doa;
    aziMeasuredDeg(iMeas) = atan2(doa(2), doa(1)) * 180 / pi;
    if aziMeasuredDeg(iMeas)<0
        aziMeasuredDeg(iMeas) = aziMeasuredDeg(iMeas) + 360;
    end
end




%% cut out direct sound

win = [ones(directSoundLength_samp-fade_samp_DS2ER-1, 1 ); ...
    cos((0:fade_samp_DS2ER)' / fade_samp_DS2ER * pi / 2).^2];

win2 = zeros(max(earlyRefsCutoff_samp)+fade_samp_ER2LR+directSoundLength_samp-fade_samp_DS2ER,numAmbiChannels,numMeasurements);
% win3 = ones(max(earlyRefsCutoff_samp)+fade_samp_ER2LR+directSoundLength_samp-fade_samp_DS2ER,numMeasurements);
win3 = ones(irLength_samp,numAmbiChannels,numMeasurements);

win2length = earlyRefsCutoff_samp+fade_samp_ER2LR+directSoundLength_samp-fade_samp_DS2ER;

for i = 1:numMeasurements
    win2(1:win2length(i),:,i) = repmat([zeros(directSoundLength_samp-fade_samp_DS2ER-1,1);...
        sin((0:fade_samp_DS2ER)' / fade_samp_DS2ER * pi / 2).^2; ...
        ones(earlyRefsCutoff_samp(i)-fade_samp_DS2ER-1, 1 ); ...
        cos((0:fade_samp_ER2LR)' / fade_samp_ER2LR * pi / 2).^2],[1,numAmbiChannels]);
    
    win3(1:win2length(i),:,i) = repmat([zeros(directSoundLength_samp-fade_samp_DS2ER-1,1);...
        zeros(earlyRefsCutoff_samp(i),1);...
        (1-cos((0:fade_samp_ER2LR)' / fade_samp_ER2LR * pi / 2).^2)],[1,numAmbiChannels]);
end

hAllDirectSound = zeros(size(hAll, 1), numAmbiChannels, numMeasurements);
% hAllEarlyRefs = zeros(size(hAll, 1), numAmbiChannels, numMeasurements);

hAllDirectSound(1:length(win),:,:) = hAll(1:length(win),1:numAmbiChannels,:) .* repmat(win,1,numAmbiChannels,numMeasurements);
% hAllEarlyRefs(1:length(win2),:,:) = hAll(1:length(win2),1:numAmbiChannels,:) .* win2;
% hAllReverb = hAll(:,1:numAmbiChannels,:) .* win3;














directSoundIdx = zeros(numMeasurements,1);
% directSoundIsolated = zeros(directSoundLength_samp, numAmbiChannels, numMeasurements);
% earlyRefsIsolated = zeros(earlyRefsCutoff_samp+fade_samp_ER2LR, numAmbiChannels, numMeasurements);
% lateRevIsolated = zeros(length(hAll(:,1,1))-earlyRefsCutoff_samp+fade_samp_DS2ER-directSoundLength_samp,numAmbiChannels, numMeasurements);

onsetThreshold = 0.5;

% for iMeas = 1:numMeasurements
%
%
%     directSoundIsolated(:, :, iMeas) = hAll(cutStart_samp:(cutStart_samp+directSoundLength_samp-1), 1:numAmbiChannels, iMeas);
%     earlyRefsIsolated(:, :, iMeas) = hAll(cutStart_samp+directSoundLength_samp-fade_samp_DS2ER:(earlyRefsCutoff_samp+directSoundLength_samp+cutStart_samp+fade_samp_ER2LR-fade_samp_DS2ER-1), 1:numAmbiChannels, iMeas);
%
%     % Cut into 3 parts:
% %     hAllDirectSound(cutStart_samp:(cutStart_samp+directSoundLength_samp-1), :, iMeas) = ...
% %         win .* directSoundIsolated(:, :, iMeas);
% %     hAllEarlyRefs(cutStart_samp+directSoundLength_samp-fade_samp_DS2ER:(earlyRefsCutoff_samp+directSoundLength_samp+fade_samp_ER2LR+cutStart_samp-fade_samp_DS2ER-1),:,iMeas) = ...
% %         win2 .* earlyRefsIsolated(:,:,iMeas);
% %     hAllReverb(cutStart_samp+directSoundLength_samp-fade_samp_DS2ER:(earlyRefsCutoff_samp+directSoundLength_samp+fade_samp_ER2LR+cutStart_samp-fade_samp_DS2ER-1), :, iMeas) = ...
% %         (1-win2) .* earlyRefsIsolated(:, :, iMeas);
% %     hAllReverb(1:cutStart_samp+directSoundLength_samp, :, iMeas) = 0;
%
%     lateRevIsolated(:,:,iMeas) = hAllReverb(earlyRefsCutoff_samp+directSoundLength_samp+fade_samp_ER2LR+cutStart_samp-fade_samp_DS2ER-fade_samp_ER2LR:end,1:numAmbiChannels,iMeas);
% end


%% Define a grid of interpolated positions

[interpolatedPositionsX_cm, interpolatedPositionsY_cm] = ...
    meshgrid(positionsInterpolatedX_cm, positionsInterpolatedY_cm);

interpolatedPositionsX_cm = interpolatedPositionsX_cm(:);
interpolatedPositionsY_cm = interpolatedPositionsY_cm(:);

numInterpolatedPoints = size(interpolatedPositionsX_cm, 1);

interpolatedPositions_cm = [interpolatedPositionsX_cm, interpolatedPositionsY_cm, ...
    receiverHeight_cm*ones(numInterpolatedPoints, 1)];

%% Spatial Filterbank / Beamforming
[~, secDirs] = getTdesign(2*NshInterp+1);
secDirs = rad2deg(secDirs);
numSecs = size(secDirs, 1);

% % Get beam weights. Cardioid, supercardioid or Max rE
% w_n = beamWeightsCardioid2Spherical(NshInterp);
% w_n = beamWeightsSupercardioid2Spherical(NshInterp);
% w_nm_re = getMaxREchannelweights(NshInterp); w_n = [w_nm_re(1); w_nm_re(2); w_nm_re(5); w_nm_re(10); w_nm_re(17)]; % use max re n weights
w_n = beamWeightsMaxEV(NshInterp); % or use the max re weights function from spherical array processing toolbox (different normalisation)

w_nm = zeros(length(hAll(1,:,1)),length(secDirs(:,1)));
for i = 1:length(secDirs(:,1))
    w_nm(:,i) = rotateAxisCoeffs(w_n, pi/2-(deg2rad(secDirs(i,2))), deg2rad(secDirs(i,1)), 'real')';
end

%% Interpolation

fax = 0:fs/nfft:fs/2;

hAllDirectSoundInterpolated = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);
hAllEarlyRefsInterpolated = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);
hAllReverbInterpolated = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);

hAllInterpolated = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);

for iInterp = 1:numInterpolatedPoints
    
    % select the closest position and keep sorted distances for interpolation
    [sortedDistances, idxSorted] = sort(vecnorm(positions - interpolatedPositions_cm(iInterp, :), 2, 2));
    idxNearest = idxSorted(1);
    
    % get start and end samples and lengths of windows for RIR portioning
    earlyRefsStart_samp = directSoundLength_samp-fade_samp_DS2ER+1;
    earlyRefsEnd_samp = earlyRefsCutoff_samp(idxNearest)+fade_samp_ER2LR+earlyRefsStart_samp-1;
    lateRevStart_samp = earlyRefsCutoff_samp(idxNearest)+earlyRefsStart_samp;
    
    nfft2 = length(earlyRefsStart_samp:earlyRefsEnd_samp);
    nfft3 = length(lateRevStart_samp:irLength_samp);
    
    earlyRefsIsolated = hAll(earlyRefsStart_samp:earlyRefsEnd_samp, 1:numAmbiChannels, idxSorted(1:2^mode_dimension)).* ...
        repmat(win2(earlyRefsStart_samp:earlyRefsEnd_samp,1:numAmbiChannels,idxNearest),[1,1,2^mode_dimension]);
    lateRevIsolated = hAll(lateRevStart_samp:end, 1:numAmbiChannels, idxSorted(1:2^mode_dimension)).* ...
        repmat(win3(lateRevStart_samp:end,1:numAmbiChannels,idxNearest),[1,1,2^mode_dimension]);
    
%% DoA interp for direct sound interpolation  
    [azi1, ele1, r] = cart2sph(directSoundDirections(idxNearest,1), ...
        directSoundDirections(idxNearest,2), ...
        directSoundDirections(idxNearest,3));
    [azi2, ele2, ~] = cart2sph(directSoundDirections(idxSorted(2),1), ...
        directSoundDirections(idxSorted(2),2), ...
        directSoundDirections(idxSorted(2),3));
    
    gains = 1 ./ (sortedDistances(1:2^mode_dimension)+eps);
    gains = gains / sum(abs(gains));
    gains2(1,1,:) = gains; % for the SH gains
    
    aziTarget = azi1 - angdiff(azi1,azi2)*gains(2);
    eleTarget = ele1 - angdiff(ele1,ele2)*gains(2);
    
    aziDif = angdiff(azi1,azi2)*gains(2);
    eleDif = angdiff(ele1,ele2)*gains(2);
    
    %% direct sound interpolation
    
    switch INTERPOLATION_MODE
        case 'rotationOnly'
            % use the nearest direct sound spectrum
            shRotMatrix = getSHrotMtx(rotz(aziCorrectionDeg), NshInterp, 'real');
            hAllDirectSoundInterpolated(:, :, iInterp) = ...
                hAllDirectSound(:, 1:numAmbiChannels, idxNearest) * shRotMatrix';
            
        case 'fixedSpectrum'
            % omits source directivity and distance by using the same direct sound
            % all the time
            idxFixed  = 4; % use this measurement for the direct sound
            shRotMatrix = getSHrotMtx(rotz(aziTarget - aziMeasuredDeg(idxFixed)), Nsh, 'real');
            hAllDirectSoundInterpolated(:, :, iMeas) = ...
                hAllDirectSound(:, 1:numAmbiChannels, idxFixed) * shRotMatrix';
            
        case 'meanSpectrum'
            % avg of FFTs
            
            if sortedDistances(1) == 0 % if it's the point of a measurement, bypass interpolation
                hAllDirectSoundInterpolated(1:directSoundLength_samp, :, iInterp) = ...
                    hAllDirectSound(1:directSoundLength_samp, 1:numAmbiChannels, idxNearest);
            else
                directSoundNearestPoints = ...
                    squeeze(hAllDirectSound(1:directSoundLength_samp,...
                    :, idxSorted(1:2^mode_dimension)));
                
                XDSNearestPoints = fft(directSoundNearestPoints);
                XDSNearestPoints_gain = XDSNearestPoints .* gains2;
                XdirectSoundNearestPointsInterp = sum(XDSNearestPoints_gain,3);
                
                hDirectInterp = real(ifft(XdirectSoundNearestPointsInterp));
                % encode to correct rotation. Remove the convert function
                % if the SRIRs are in N3D normalisation
                hDirectEncoded = convert_N3D_SN3D(rotateHOA_N3D(convert_N3D_SN3D(hDirectInterp,'sn2n'),aziDif,eleDif,0),'n2sn');
                
                % normalise RMS of direct sound
                rmsNormValue = sum(rms(squeeze(directSoundNearestPoints(:,1,:)))*gains);
                hDirectEncoded = hDirectEncoded / rms(hDirectEncoded(:, 1)) * rmsNormValue;
                hAllDirectSoundInterpolated(1:directSoundLength_samp, :, iInterp) = ...
                    hDirectEncoded;
            end
            
        case 'minPhase'
            % interpolate direct sound spectrum of the nearest 4 points
            % (if 2d) or 2 points (if 1d)
            
            directSoundNearestPoints = ...
                squeeze(hAllDirectSound(1:directSoundLength_samp, 1, idxSorted(1:2^mode_dimension)));
            XdirectSoundNearestPoints = fft(directSoundNearestPoints, nfft);
            XdirectSoundNearestPointsSmooth = ...
                smoothSpectrum(abs(XdirectSoundNearestPoints(1:nfft/2+1, :)), fax(:), 3);
            XdirectSoundNearestPointsInterp = XdirectSoundNearestPointsSmooth * gains;
            hMinPhase = designMinPhase(XdirectSoundNearestPointsInterp);
            
            % encode to correct rotation. Remove the convert function
            % if the SRIRs are in N3D normalisation
            hDirectEncoded = hMinPhase(1:directSoundLength_samp) * ...
                convert_N3D_SN3D(evalSH(NshInterp, [aziTarget, eleTarget]),'n2sn');
            
            % normalise RMS of direct sound
            rmsNormValue = sum(rms(directSoundNearestPoints)*gains);
            hDirectEncoded = hDirectEncoded / rms(hDirectEncoded(:, 1)) * rmsNormValue;
            hAllDirectSoundInterpolated(1:directSoundLength_samp, :, iInterp) = ...
                hDirectEncoded;
    end
    
    
    %% interpolate early reflections
    % interpolate direct sound spectrum of the nearest 4 points if 2d
    % and 2 points if 1d
    if sortedDistances(1) == 0 % if it's an exact point, don't bother interpolating
        hAllEarlyRefsInterpolated(earlyRefsStart_samp:earlyRefsEnd_samp,:,iInterp) = ...
        earlyRefsIsolated(:,1:numAmbiChannels,1);
%             hAllEarlyRefs(:, 1:numAmbiChannels, idxNearest);
    else
        % take fft of nearest measurements, weight by gains
%         earlyNearestPoints = earlyRefsIsolated;
        XearlyRefsNearestPoints = fft(earlyRefsIsolated);
        XearlyRefsNearestPoints_gain = XearlyRefsNearestPoints .* gains2;
        
        % then get the weighted nearest measurements in sectors
        rir_secs = zeros(length(XearlyRefsNearestPoints_gain(:,1,1)),numSecs,length(XearlyRefsNearestPoints_gain(1,1,:)));
        for i = 1:2^mode_dimension
            rir_secs(:,:,i) = XearlyRefsNearestPoints_gain(:,:,i) * w_nm;
        end
        
        rir_secs_interp = zeros(length(XearlyRefsNearestPoints_gain(:,1,1)),numSecs);
        
        for secIdx = 1:numSecs % for each sector:
            rir_channelNearests = rir_secs(:, secIdx,:);
            
            XearlyRefsNearestPointsInterp = sum(rir_channelNearests,3);
            
            % over frequency bands, obtain the RMS magnitude difference of the
            % interpolated and nearest IRs (as comb filtering etc produces
            % reductions in magnitude). currently lin freq bands because lazy
            % (should probably be ERB)
            n_bands = nfft2/20; % must be integer
            correctRMS = zeros(n_bands,1);
            currentRMS = zeros(n_bands,1);
            for j = 1:n_bands
                % calculate rms of freq bands
                correctRMS(j,:) = sum(rms(rir_channelNearests(round((j-1)*nfft2/n_bands+1):round(j*nfft2/n_bands),1,:)),3);
                currentRMS(j,:) = rms(XearlyRefsNearestPointsInterp(round((j-1)*nfft2/n_bands+1):round(j*nfft2/n_bands),1));
            end
            rmsDiff = correctRMS./currentRMS;
            
            % create array of RMS equalisation gains
            rmsGains = zeros(length(XearlyRefsNearestPointsInterp),1);
            for j = 1:n_bands
                if j == n_bands % last band
                    rmsGains((j-1)*nfft2/n_bands+1:j*nfft2/n_bands,1) = linspace(rmsDiff(j,1),rmsDiff(j,1),length(round((j-1)*nfft2/n_bands+1):round(j*nfft2/n_bands)));
                else
                    rmsGains((j-1)*nfft2/n_bands+1:j*nfft2/n_bands,1) = linspace(rmsDiff(j,1),rmsDiff(j+1,1),length(round((j-1)*nfft2/n_bands+1):round(j*nfft2/n_bands)));
                end
            end
            % apply gains
            HEarlyMinPhaseSec = XearlyRefsNearestPointsInterp.*rmsGains;
            rir_secs_interp(:,secIdx) = HEarlyMinPhaseSec;
        end
        
        % back into SH domain (still in Freq domain though)
        HEarlyMinPhase = rir_secs_interp  * w_nm'; % my own sector stuff
        
        % needs some normalisation due to the beamWeightsMaxEV stuff
        w_n_gain = w_n.^2 * length(secDirs(:,1)) ./ ((0:1:NshInterp) * 2 + 1)';
        for i = 1:length(HEarlyMinPhase(1,:))
            n = ceil(sqrt(i)-1);
            HEarlyMinPhase(:,i) = HEarlyMinPhase(:,i) / w_n_gain(n+1);
        end
        
        % return to time-domain and window
        hEarlyMinPhase = real(ifft(HEarlyMinPhase)) ;%.* win2;
        
        % normalise interpolated early reflections -- for each SH channel
        rmsNormValueER = sum(rms(earlyRefsIsolated) .* repmat(gains2,1,numAmbiChannels),3);
        hEarlyMinPhase = hEarlyMinPhase ./ rms(hEarlyMinPhase) .* rmsNormValueER;
        
        hAllEarlyRefsInterpolated(earlyRefsStart_samp:earlyRefsEnd_samp, :, iInterp) = ...
            hEarlyMinPhase;
    end
    
    %% Interpolate late reverb
    if sortedDistances(1) == 0 % if it's an exact point, don't bother interpolating
        hAllReverbInterpolated(lateRevStart_samp:end,:,iInterp) = ...
            lateRevIsolated(:, 1:numAmbiChannels, 1);
%             hAllReverb(:, 1:numAmbiChannels, idxNearest);
    else
%         lateNearestPoints = ...
%             squeeze(lateRevIsolated(:, :, idxSorted(1:2^mode_dimension)));
        XlateNearestPoints = fft(lateRevIsolated);
        
        XlateNearestPoints_gain = XlateNearestPoints .* gains2;
        XlateNearestPointsInterp = sum(XlateNearestPoints_gain,3);
        
        % over frequency bands, obtain the RMS magnitude difference of the
        % interpolated and nearest IRs (as comb filtering etc produces
        % reductions in magnitude). currently lin freq bands because lazy
        % (should probably be ERB)
        n_bands = nfft3/20; % must be integer
        correctRMS_late = zeros(n_bands,1);
        currentRMS_late = zeros(n_bands,1);
        for j = 1:n_bands
            % calculate rms of freq bands
            correctRMS_late(j) = sum(rms((XlateNearestPoints_gain(round((j-1)*(nfft3/n_bands)+1):round((j)*(nfft3/n_bands)),1,:))),3);
            currentRMS_late(j) = rms((XlateNearestPointsInterp(round((j-1)*(nfft3/n_bands)+1):round((j)*(nfft3/n_bands)),1)));
        end
        rmsDiff_late = correctRMS_late./currentRMS_late;
        
        % create matrix of RMS equalisation gains
        rmsGains_late = zeros(length(XlateNearestPointsInterp),1);
        for j = 1:n_bands
            if j == n_bands % last band
                rmsGains_late(((j-1)*(nfft3/n_bands)+1):((j)*(nfft3/n_bands)),1) = linspace(rmsDiff_late(j,1),rmsDiff_late(j,1),length(round((j-1)*(nfft3/n_bands)+1):round((j)*(nfft3/n_bands))));
            else
                rmsGains_late(((j-1)*(nfft3/n_bands)+1):((j)*(nfft3/n_bands)),1) = linspace(rmsDiff_late(j,1),rmsDiff_late(j+1,1),length(round((j-1)*(nfft3/n_bands)+1):round((j)*(nfft3/n_bands))));
            end
        end
        % apply gains, return to time domain and windowing at start and end of ERs
        hLateMinPhase = real(ifft(XlateNearestPointsInterp.*repmat(rmsGains_late,[1,numAmbiChannels]))); % .* repmat(win2,[1,numAmbiChannels]);
        
        % normalise late reverb for each SH channel
        rmsNormValueLR = sum(rms(lateRevIsolated) .* repmat(gains2,1,numAmbiChannels),3);
        hLateMinPhase = hLateMinPhase ./ rms(hLateMinPhase) .* rmsNormValueLR;
        
        hAllReverbInterpolated(lateRevStart_samp:end, :, iInterp) = ...
            hLateMinPhase;
    end
    
    %% Construct the final IRs
    
    %     % % Bypass direct sound interpolation (debug)
    %     hAllDirectSoundInterpolated(:,:,iInterp) = ...
    %         hAllDirectSound(:, 1:numAmbiChannels, idxNearest);
    %
    %     % % Bypass ER interpolation (debug)
    %     hAllEarlyRefsInterpolated(:,:,iInterp) = ...
    %         hAllEarlyRefs(:, 1:numAmbiChannels, idxNearest);
    %
    %     % % Bypass LR interpolation (debug)
    %     hAllReverbInterpolated(:,:,iInterp) = ...
    %         hAllReverb(:, 1:numAmbiChannels, idxNearest);
    
    % add three parts together
    hAllInterpolated(:, :, iInterp) = hAllDirectSoundInterpolated(:, :, iInterp) * 10^(directMod_db/20)+ ...
        hAllEarlyRefsInterpolated(:, :, iInterp) * 10^(roomMod_db/20) +...
        hAllReverbInterpolated(:, :, iInterp) * 10^(roomMod_db/20);
    
end

%%
interpType = 'sectors, no time windows';
run analyse_interpolated_SRIRs.m
% figure;plot(indexMin)
%

% figure;plot(r1);hold on; plot(rmsNormValLR_);plot(rmsNormValLR_2);title('lr rms after mix / norm value / before mix');

% figure;hold on;%for i = 1:7:7
% plot(hAll(1:2000,1,1));
% plot(hAll(1:2000,1,7));
% %end
%% Save to sofa
% %{
SOFAstart()

Obj = sofa;
Obj.Data.IR = permute(hAllInterpolated, [3 2 1]);

% Update dimensions
Obj=SOFAupdateDimensions(Obj);

% Fill with attributes
Obj.GLOBAL_ListenerShortName = 'EM';
Obj.GLOBAL_History = 'created on 07.03.2022';
Obj.GLOBAL_DatabaseName = 'none';
Obj.GLOBAL_ApplicationName = 'SOFA API';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
Obj.GLOBAL_Organization = 'Aalto Acoustics Lab';
Obj.GLOBAL_AuthorContact = 'thomas.mckenzie@aalto.fi';
Obj.GLOBAL_Comment = ' ';
Obj.GLOABL_Title =  'Responses on a line, interpolated';

% save the SOFA file
switch INTERPOLATION_MODE
    case 'meanSpectrum'
        SOFAfn = fullfile(['srirInterp_ms_sectors_winlen_',num2str(measReduction),'.sofa']);
    case 'minPhase'
        SOFAfn = fullfile(['srirInterp_mp_',num2str(measReduction),'.sofa']);
end

disp(['Saving:  ' SOFAfn]);
Obj = SOFAsave(SOFAfn, Obj, 1);
%}
% end