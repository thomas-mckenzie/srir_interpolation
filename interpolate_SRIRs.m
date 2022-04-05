function [srirs_interp,pos_interp] = interpolate_SRIRs(srirs_input,pos_input,resolution_new,fs,interp_modeDS,length_sampDS,fade_sampDS2ER,fade_sampER2LR,N_interp)
% Perceptual interpolation of SRIRs
%   srirs_input - input format [samples, sh_channels, num_measurements]

if nargin<9; N_interp = sqrt(size(srirs_input,2))-1;    end
if nargin<8; fade_sampER2LR = fs/100;                   end
if nargin<7; fade_sampDS2ER = 20;                       end
if nargin<6; length_sampDS = 200;                       end
if nargin<5; interp_modeDS = 'minPhase';                end

%% init

irLength_samp = size(srirs_input,1);
numMeasurements = size(srirs_input,3);
numAmbiChannels = size(srirs_input,2);

nfftDS = 2^nextpow2(length_sampDS);

%% Interpolated positions

% work out whether interpolation should be 1D, 2D or 3D
isPlaneUsed=zeros(size(pos_input,2),1);
for i = 1:size(pos_input,2)
    if all(pos_input(:,i) == pos_input(1,i))
        isPlaneUsed(i) = 0;
    else
        isPlaneUsed(i) = 1;
    end
end
mode_dimension = sum(isPlaneUsed); % 1 = 1D (planar), 2 = 2D, 3 = 3D

% calculate the interp positions
if mode_dimension == 1 % 1D - assumes interpolating on X axis
    positionsInterpolatedX_cm = min(pos_input(:,1)):resolution_new:max(pos_input(:,1));
    positionsInterpolatedY_cm = pos_input(1,2);
    positionsInterpolatedZ_cm = pos_input(1,3);
elseif mode_dimension == 2 % 2D - assumes interpolating on X and Y axes
    positionsInterpolatedX_cm = min(pos_input(:,1)):resolution_new:max(pos_input(:,1));
    positionsInterpolatedY_cm = min(pos_input(:,2)):resolution_new:max(pos_input(:,2));
    positionsInterpolatedZ_cm = pos_input(1,3);   
else % 3D
    positionsInterpolatedX_cm = min(pos_input(:,1)):resolution_new:max(pos_input(:,1));
    positionsInterpolatedY_cm = min(pos_input(:,2)):resolution_new:max(pos_input(:,2));
    positionsInterpolatedZ_cm = min(pos_input(:,3)):resolution_new:max(pos_input(:,3));
end

% Define a grid of interpolated positions
[interpolatedPositionsX_cm, interpolatedPositionsY_cm, interpolatedPositionsZ_cm] = ...
    meshgrid(positionsInterpolatedX_cm, positionsInterpolatedY_cm, positionsInterpolatedZ_cm);

interpolatedPositionsX_cm = interpolatedPositionsX_cm(:);
interpolatedPositionsY_cm = interpolatedPositionsY_cm(:);
interpolatedPositionsZ_cm = interpolatedPositionsZ_cm(:);

pos_interp = [interpolatedPositionsX_cm, interpolatedPositionsY_cm, ...
    interpolatedPositionsZ_cm];
numInterpolatedPoints = length(pos_interp);

%% Calculate the cutoff between early refs and late reverb
% based on when the energy decay curve amplitude is lower than a target
% in Hikada a similar method uses 1/exp(1) but this was found to be too
% early.

% target = 1 / exp(1);
ER_target = 0.1; % empirically chosen value.
earlyRefsCutoff_samp = zeros(1,numMeasurements);
for i = 1:numMeasurements
    decay = rir2decay(srirs_input(:,1,i), fs, 1000, 1, 1, 1);
    ind_at_t_L = find(decay<ER_target,1);
    
    earlyRefsCutoff_samp(i) = ceil(ind_at_t_L/1000) * 1000; % nearest 1000 (empirically chosen)
end

%% initial DOA (used for interpolating direct sound)

m = SRIR_Measurement();m.setShMicArray(4, 0.032);m.fs = fs; % for eigenmike
directSoundDirections = zeros(numMeasurements, 3);
aziMeasuredDeg = zeros(numMeasurements, 1);
for iMeas = 1:numMeasurements
    m.srir = srirs_input(1:length_sampDS, :, iMeas);
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

%% cut out direct sound and make windows for early refs and late reverberation
% windows for ER and LR can vary with the ER window size (earlyRefsCutoff_samp)

winDS = [ones(length_sampDS-fade_sampDS2ER-1, 1); cos((0:fade_sampDS2ER)' / fade_sampDS2ER * pi / 2).^2];
winER = zeros(max(earlyRefsCutoff_samp)+fade_sampER2LR+length_sampDS-fade_sampDS2ER,numAmbiChannels,numMeasurements);
winLR = ones(irLength_samp,numAmbiChannels,numMeasurements);

winERlength = earlyRefsCutoff_samp+fade_sampER2LR+length_sampDS-fade_sampDS2ER;
for i = 1:numMeasurements
    winER(1:winERlength(i),:,i) = repmat([zeros(length_sampDS-fade_sampDS2ER-1,1);...
        sin((0:fade_sampDS2ER)' / fade_sampDS2ER * pi / 2).^2; ...
        ones(earlyRefsCutoff_samp(i)-fade_sampDS2ER-1, 1 ); ...
        cos((0:fade_sampER2LR)' / fade_sampER2LR * pi / 2).^2],[1,numAmbiChannels]);
    
    winLR(1:winERlength(i),:,i) = repmat([zeros(length_sampDS-fade_sampDS2ER-1,1);...
        zeros(earlyRefsCutoff_samp(i),1);...
        (1-cos((0:fade_sampER2LR)' / fade_sampER2LR * pi / 2).^2)],[1,numAmbiChannels]);
end

hAllDirectSound = zeros(size(srirs_input, 1), numAmbiChannels, numMeasurements);
hAllDirectSound(1:length(winDS),:,:) = srirs_input(1:length(winDS),1:numAmbiChannels,:) .* repmat(winDS,1,numAmbiChannels,numMeasurements);

%% Spherical Filter Bank / Beamforming
[~, secDirs] = getTdesign(2*N_interp+1);
secDirs = rad2deg(secDirs);
numSecs = size(secDirs, 1);

% % Get beam weights. Cardioid, supercardioid or Max rE
% w_n = beamWeightsCardioid2Spherical(NshInterp);
% w_n = beamWeightsSupercardioid2Spherical(NshInterp);
% w_nm_re = getMaxREchannelweights(NshInterp); w_n = [w_nm_re(1); w_nm_re(2); w_nm_re(5); w_nm_re(10); w_nm_re(17)]; % use max re n weights
w_n = beamWeightsMaxEV(N_interp); % or use the max re weights function from spherical array processing toolbox (different normalisation)

w_nm = zeros(size(srirs_input, 2),size(secDirs,1));
for i = 1:length(secDirs(:,1))
    w_nm(:,i) = rotateAxisCoeffs(w_n, pi/2-(deg2rad(secDirs(i,2))), deg2rad(secDirs(i,1)), 'real')';
end

%% Interpolation

fax = 0:fs/nfftDS:fs/2;
hAllDirectSoundInterpolated = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);
hAllEarlyRefsInterpolated = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);
hAllReverbInterpolated = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);
srirs_interp = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);

for iInterp = 1:numInterpolatedPoints
    % calculate closest measurements and distances for interpolation
    [sortedDistances, idxSorted] = sort(vecnorm(pos_input - pos_interp(iInterp, :), 2, 2));
    idxNearest = idxSorted(1); % index of closest measurement
    
    % get start and end samples and lengths of windows for RIR portioning
    earlyRefsStart_samp = length_sampDS-fade_sampDS2ER+1;
    earlyRefsEnd_samp = earlyRefsCutoff_samp(idxNearest)+fade_sampER2LR+earlyRefsStart_samp-1;
    lateRevStart_samp = earlyRefsCutoff_samp(idxNearest)+earlyRefsStart_samp;
    
    % fft window sizes
    nfftER = length(earlyRefsStart_samp:earlyRefsEnd_samp);
    nfftLR = length(lateRevStart_samp:irLength_samp);
    
    % isolate early and late reverb, and window
    earlyRefsIsolated = srirs_input(earlyRefsStart_samp:earlyRefsEnd_samp, 1:numAmbiChannels, idxSorted(1:2^mode_dimension)).* ...
        repmat(winER(earlyRefsStart_samp:earlyRefsEnd_samp,1:numAmbiChannels,idxNearest),[1,1,2^mode_dimension]);
    lateRevIsolated = srirs_input(lateRevStart_samp:end, 1:numAmbiChannels, idxSorted(1:2^mode_dimension)).* ...
        repmat(winLR(lateRevStart_samp:end,1:numAmbiChannels,idxNearest),[1,1,2^mode_dimension]);
    
    %% interpolate direct sound
    % DoA interp for direct sound interpolation
    [azi1, ele1, r] = cart2sph(directSoundDirections(idxNearest,1), ...
        directSoundDirections(idxNearest,2), ...
        directSoundDirections(idxNearest,3));
    [azi2, ele2, ~] = cart2sph(directSoundDirections(idxSorted(2),1), ...
        directSoundDirections(idxSorted(2),2), ...
        directSoundDirections(idxSorted(2),3));
    
    gains = 1 ./ (sortedDistances(1:2^mode_dimension)+eps);
    gains = gains / sum(abs(gains));
    gainsSH(1,1,:) = gains; % for the SH channel gains (ER and LR)
    
    aziTarget = azi1 - angdiff(azi1,azi2)*gains(2);
    eleTarget = ele1 - angdiff(ele1,ele2)*gains(2);
    aziDif = angdiff(azi1,azi2)*gains(2);
    eleDif = angdiff(ele1,ele2)*gains(2);
    
    switch interp_modeDS
        case 'rotationOnly'
            % use the nearest direct sound spectrum
            shRotMatrix = getSHrotMtx(rotz(aziCorrectionDeg), N_interp, 'real');
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
                hAllDirectSoundInterpolated(1:length_sampDS, :, iInterp) = ...
                    hAllDirectSound(1:length_sampDS, 1:numAmbiChannels, idxNearest);
            else
                directSoundNearestPoints = ...
                    squeeze(hAllDirectSound(1:length_sampDS,...
                    :, idxSorted(1:2^mode_dimension)));
                
                XDSNearestPoints = fft(directSoundNearestPoints);
                XDSNearestPoints_gain = XDSNearestPoints .* gainsSH;
                XdirectSoundNearestPointsInterp = sum(XDSNearestPoints_gain,3);
                
                hDirectInterp = real(ifft(XdirectSoundNearestPointsInterp));
                % encode to correct rotation. Remove the convert function
                % if the SRIRs are already in N3D normalisation
                hDirectEncoded = convert_N3D_SN3D(rotateHOA_N3D(convert_N3D_SN3D(hDirectInterp,'sn2n'),aziDif,eleDif,0),'n2sn');
                
                % normalise RMS of direct sound
                rmsNormValue = sum(rms(squeeze(directSoundNearestPoints(:,1,:))) * gains);
                hDirectEncoded = hDirectEncoded / rms(hDirectEncoded(:, 1)) * rmsNormValue;
                hAllDirectSoundInterpolated(1:length_sampDS, :, iInterp) = ...
                    hDirectEncoded;
            end
            
        case 'minPhase'
            % interpolate direct sound spectrum of the nearest 4 points
            % (if 2d) or 2 points (if 1d)
            
            directSoundNearestPoints = ...
                squeeze(hAllDirectSound(1:length_sampDS, 1, idxSorted(1:2^mode_dimension)));
            XdirectSoundNearestPoints = fft(directSoundNearestPoints, nfftDS);
            XdirectSoundNearestPointsSmooth = ...
                smoothSpectrum(abs(XdirectSoundNearestPoints(1:nfftDS/2+1, :)), fax(:), 3);
            XdirectSoundNearestPointsInterp = XdirectSoundNearestPointsSmooth * gains;
            hMinPhase = designMinPhase(XdirectSoundNearestPointsInterp);
            
            % encode to correct rotation. Remove the convert function
            % if the SRIRs are in N3D normalisation
            hDirectEncoded = hMinPhase(1:length_sampDS) * ...
                convert_N3D_SN3D(evalSH(N_interp, [aziTarget, eleTarget]),'n2sn');
            
            % normalise RMS of direct sound
            rmsNormValue = sum(rms(squeeze(directSoundNearestPoints(:,1,:))) * gains);
            hDirectEncoded = hDirectEncoded / rms(hDirectEncoded(:, 1)) * rmsNormValue;
            hAllDirectSoundInterpolated(1:length_sampDS, :, iInterp) = ...
                hDirectEncoded;
    end
    
    
    %% interpolate early reflections
    if sortedDistances(1) == 0 % if it's an exact point, don't bother interpolating
        hAllEarlyRefsInterpolated(earlyRefsStart_samp:earlyRefsEnd_samp,:,iInterp) = ...
            earlyRefsIsolated(:,1:numAmbiChannels,1);
    else
        % take fft of nearest measurements, weight by gains
        XearlyRefsNearestPoints = fft(earlyRefsIsolated);
        XearlyRefsNearestPoints_gain = XearlyRefsNearestPoints .* gainsSH;
        
        % then get the weighted nearest measurements in sectors
        XearlyRefsNearestPoints_gain_secs = zeros(length(XearlyRefsNearestPoints_gain(:,1,1)),numSecs,length(XearlyRefsNearestPoints_gain(1,1,:)));
        for i = 1:2^mode_dimension
            XearlyRefsNearestPoints_gain_secs(:,:,i) = XearlyRefsNearestPoints_gain(:,:,i) * w_nm;
        end
        
        HEarlyInterp_secs = zeros(length(XearlyRefsNearestPoints_gain(:,1,1)),numSecs);
        for secIdx = 1:numSecs % for each sector:
            XearlyRefsNearestPoints_gain_sector = XearlyRefsNearestPoints_gain_secs(:, secIdx,:);
            XearlyRefsNearestPointsInterp_sector = sum(XearlyRefsNearestPoints_gain_sector,3);
            
            % over frequency bands, obtain the RMS magnitude difference of the
            % interpolated and nearest IRs (as comb filtering etc produces reductions in magnitude). 
            n_bandsER = nfftER/20; % must be integer
            rmsER_target = zeros(n_bandsER,1);
            rmsER_current = zeros(n_bandsER,1);
            for j = 1:n_bandsER
                % calculate rms of freq bands
                rmsER_target(j,:) = sum(rms(XearlyRefsNearestPoints_gain_sector(round((j-1)*nfftER/n_bandsER+1):round(j*nfftER/n_bandsER),1,:)),3);
                rmsER_current(j,:) = rms(XearlyRefsNearestPointsInterp_sector(round((j-1)*nfftER/n_bandsER+1):round(j*nfftER/n_bandsER),1));
            end
            rmsDiffER = rmsER_target./rmsER_current;
            
            % create array of RMS equalisation gains
            rmsGainsER = zeros(length(XearlyRefsNearestPointsInterp_sector),1);
            for j = 1:n_bandsER
                if j == n_bandsER % last band
                    rmsGainsER((j-1)*nfftER/n_bandsER+1:j*nfftER/n_bandsER,1) = ...
                        linspace(rmsDiffER(j,1),rmsDiffER(j,1),length(round((j-1)*nfftER/n_bandsER+1):round(j*nfftER/n_bandsER)));
                else
                    rmsGainsER((j-1)*nfftER/n_bandsER+1:j*nfftER/n_bandsER,1) = ...
                        linspace(rmsDiffER(j,1),rmsDiffER(j+1,1),length(round((j-1)*nfftER/n_bandsER+1):round(j*nfftER/n_bandsER)));
                end
            end
            % apply gains
            HEarlyInterp_sector = XearlyRefsNearestPointsInterp_sector.*rmsGainsER;
            HEarlyInterp_secs(:,secIdx) = HEarlyInterp_sector;
        end
        
        % back into SH domain (still in Freq domain though)
        HEarlyInterp = HEarlyInterp_secs  * w_nm'; % my own sector stuff
        
        % needs some normalisation due to the beamWeightsMaxEV stuff
        w_n_gain = w_n.^2 * length(secDirs(:,1)) ./ ((0:1:N_interp) * 2 + 1)';
        for i = 1:length(HEarlyInterp(1,:))
            n = ceil(sqrt(i)-1);
            HEarlyInterp(:,i) = HEarlyInterp(:,i) / w_n_gain(n+1);
        end
        
        % return to time-domain and window
        hEarlyInterp = real(ifft(HEarlyInterp)) ;%.* win2;
        
        % normalise interpolated early reflections -- for each SH channel
        rmsNormValueER = sum(rms(earlyRefsIsolated) .* repmat(gainsSH,1,numAmbiChannels),3);
        hEarlyInterp = hEarlyInterp ./ rms(hEarlyInterp) .* rmsNormValueER;
        
        hAllEarlyRefsInterpolated(earlyRefsStart_samp:earlyRefsEnd_samp, :, iInterp) = ...
            hEarlyInterp;
    end
    
    %% Interpolate late reverb
    if sortedDistances(1) == 0 % if it's an exact point, don't bother interpolating
        hAllReverbInterpolated(lateRevStart_samp:end,:,iInterp) = ...
            lateRevIsolated(:, 1:numAmbiChannels, 1);
    else
        XlateNearestPoints = fft(lateRevIsolated);
        
        XlateNearestPoints_gain = XlateNearestPoints .* gainsSH;
        XlateNearestPointsInterp = sum(XlateNearestPoints_gain,3);
        
        % over frequency bands, obtain the RMS magnitude difference of the
        % interpolated and nearest IRs (as comb filtering etc produces reductions in magnitude). 
        n_bandsLR = nfftLR/20; % must be integer
        rmsLR_target = zeros(n_bandsLR,1);
        rmsLR_current = zeros(n_bandsLR,1);
        for j = 1:n_bandsLR
            % calculate rms of freq bands
            rmsLR_target(j) = sum(rms((XlateNearestPoints_gain(round((j-1)*(nfftLR/n_bandsLR)+1):round((j)*(nfftLR/n_bandsLR)),1,:))),3);
            rmsLR_current(j) = rms((XlateNearestPointsInterp(round((j-1)*(nfftLR/n_bandsLR)+1):round((j)*(nfftLR/n_bandsLR)),1)));
        end
        rmsDiffLR = rmsLR_target./rmsLR_current;
        
        % create matrix of RMS equalisation gains
        rmsGainsLR = zeros(length(XlateNearestPointsInterp),1);
        for j = 1:n_bandsLR
            if j == n_bandsLR % last band
                rmsGainsLR(((j-1)*(nfftLR/n_bandsLR)+1):((j)*(nfftLR/n_bandsLR)),1) = ...
                    linspace(rmsDiffLR(j,1),rmsDiffLR(j,1),length(round((j-1)*(nfftLR/n_bandsLR)+1):round((j)*(nfftLR/n_bandsLR))));
            else
                rmsGainsLR(((j-1)*(nfftLR/n_bandsLR)+1):((j)*(nfftLR/n_bandsLR)),1) = ...
                    linspace(rmsDiffLR(j,1),rmsDiffLR(j+1,1),length(round((j-1)*(nfftLR/n_bandsLR)+1):round((j)*(nfftLR/n_bandsLR))));
            end
        end
        % apply gains, return to time domain and windowing at start and end of ERs
        hLateInterp = real(ifft(XlateNearestPointsInterp.*repmat(rmsGainsLR,[1,numAmbiChannels]))); % .* repmat(win2,[1,numAmbiChannels]);
        
        % normalise late reverb for each SH channel
        rmsNormValueLR = sum(rms(lateRevIsolated) .* repmat(gainsSH,1,numAmbiChannels),3);
        hLateInterp = hLateInterp ./ rms(hLateInterp) .* rmsNormValueLR;
        
        hAllReverbInterpolated(lateRevStart_samp:end, :, iInterp) = ...
            hLateInterp;
    end
    
    %% Construct the final IRs
    % add three parts together
    srirs_interp(:, :, iInterp) = hAllDirectSoundInterpolated(:, :, iInterp) + ...
        hAllEarlyRefsInterpolated(:, :, iInterp) +...
        hAllReverbInterpolated(:, :, iInterp);
    
end
end