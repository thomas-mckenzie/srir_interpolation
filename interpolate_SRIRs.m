function [srirs_interp,pos_interp] = interpolate_SRIRs(srirs_input,pos_input,resolution_new,fs,interp_modeDS,length_sampDS,fade_sampDS2ER,fade_sampER2LR,N_interp)
% Perceptual interpolation of SRIRs
% interpolate_SRIRs()
% srirs_input - input format [samples, sh_channels, num_measurements]
% 
% The user is directed to the ICA 2022 paper for details:
% ﻿McKenzie, T., Meyer-Kahlen, N., Daugintis, R., McCormack, L., Schlecht, S. 
% J., & Pulkki, V. (2022). Perceptual interpolation and rendering of coupled 
% room spatial room impulse responses. International Congress on Acoustics, 
% Korea. 
% 
% Thomas McKenzie, 2022. thomas.mckenzie@aalto.fi / tom.mckenzie07@gmail.com

if nargin<9; N_interp = sqrt(size(srirs_input,2))-1;    end
if nargin<8; fade_sampER2LR = fs/100;                   end
if nargin<7; fade_sampDS2ER = 20;                       end
if nargin<6; length_sampDS = 200;                       end
if nargin<5; interp_modeDS = 'minPhase';                end

%% init

irLength_samp = size(srirs_input,1);
numMeasurements = size(srirs_input,3);
numAmbiChannels = (N_interp+1)^2;

nfftDS = 2^nextpow2(length_sampDS);

% ERB parameters (Glasberg and Moore)
erb_Q = 9.26449;
erb_min = 24.7;
erb_low_freq = 10;
n_bandsER = 48;
n_bandsLR = 48;

%% Interpolated positions

% work out whether interpolation should be 1D, 2D or 3D
isPlaneUsed=zeros(size(pos_input,2),1);
for i = 1:size(pos_input,2)
    if all(pos_input(:,i) == pos_input(1,i)); isPlaneUsed(i) = 0;
    else; isPlaneUsed(i) = 1; end
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

ER_target = 0.1; % empirically chosen value. could be 1 / exp(1);
earlyRefsCutoff_samp = zeros(1,numMeasurements);
[B_decay, A_decay] = octdsgn(1000, fs, 3);
for i = 1:numMeasurements
    % filter and schroeder integration to get EDC. W channel, 1kHz freq band
    srir_input_filt = filter(B_decay,A_decay,srirs_input(:,1,i));
    int_sch = 1/length(srir_input_filt) * cumtrapz(flip(srir_input_filt/max(abs(srir_input_filt))).^2);
    edc = flip(int_sch(2:end));
    edc = edc/max(edc);
    
    ind_at_t_L = find(edc<ER_target,1);
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

hAllDS = zeros(size(srirs_input, 1), numAmbiChannels, numMeasurements);
hAllDS(1:length(winDS),:,:) = srirs_input(1:length(winDS),1:numAmbiChannels,:) .* repmat(winDS,1,numAmbiChannels,numMeasurements);

%% Spherical Filter Bank / Beamforming
[~, secDirs] = getTdesign(2*N_interp+1);
secDirs = rad2deg(secDirs);
numSecs = size(secDirs, 1);

% Get beam weights. Using max energy vector weights from spherical array processing toolbox
w_n = beamWeightsMaxEV(N_interp);

w_nm = zeros(numAmbiChannels,numSecs);
for i = 1:length(secDirs(:,1))
    w_nm(:,i) = rotateAxisCoeffs(w_n, pi/2-(deg2rad(secDirs(i,2))), deg2rad(secDirs(i,1)), 'real')';
end

%% Interpolation

fax = 0:fs/nfftDS:fs/2;
hAllDS_interp = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);
hAllER_interp = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);
hAllLR_interp = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);
srirs_interp = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);

for iInterp = 1:numInterpolatedPoints
    % calculate closest measurements and distances for interpolation
    [sortedDistances, idxSorted] = sort(vecnorm(pos_input - pos_interp(iInterp, :), 2, 2));
    idxNearest = idxSorted(1); % index of closest measurement
    
    gains = 1 ./ (sortedDistances(1:2^mode_dimension)+eps);
    gains = gains / sum(abs(gains));
    gainsSH(1,1,:) = gains; % for the SH channel gains (ER and LR)
    
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
    [azi1,ele1,~] = cart2sph(directSoundDirections(idxNearest,1), ...
        directSoundDirections(idxNearest,2), ...
        directSoundDirections(idxNearest,3));
    [azi2,ele2,~] = cart2sph(directSoundDirections(idxSorted(2),1), ...
        directSoundDirections(idxSorted(2),2), ...
        directSoundDirections(idxSorted(2),3));
    
    aziTarget = azi1 - angdiff(azi1,azi2)*gains(2);
    eleTarget = ele1 - angdiff(ele1,ele2)*gains(2);
    aziDif = angdiff(azi1,azi2)*gains(2);
    eleDif = angdiff(ele1,ele2)*gains(2);
    
    switch interp_modeDS
        case 'rotationOnly'
            % use the nearest direct sound spectrum
            shRotMatrix = getSHrotMtx(rotz(aziCorrectionDeg), N_interp, 'real');
            hAllDS_interp(:, :, iInterp) = ...
                hAllDS(:, 1:numAmbiChannels, idxNearest) * shRotMatrix';
            
        case 'fixedSpectrum'
            % omits source directivity and distance by using the same direct sound
            % all the time
            idxFixed  = 4; % use this measurement for the direct sound
            shRotMatrix = getSHrotMtx(rotz(aziTarget - aziMeasuredDeg(idxFixed)), Nsh, 'real');
            hAllDS_interp(:,:,iMeas) = ...
                hAllDS(:, 1:numAmbiChannels, idxFixed) * shRotMatrix';
            
        case 'meanSpectrum'
            % avg of FFTs
            if sortedDistances(1) == 0 % if it's the point of a measurement, bypass interpolation
                hAllDS_interp(1:length_sampDS, :, iInterp) = ...
                    hAllDS(1:length_sampDS, 1:numAmbiChannels, idxNearest);
            else
                nearestMeasDS = ...
                    squeeze(hAllDS(1:length_sampDS, :, idxSorted(1:2^mode_dimension)));
                
                XDSnearestMeas = fft(nearestMeasDS);
                XDSnearestMeas_gain = XDSnearestMeas .* gainsSH;
                XDSnearestMeas_interp = sum(XDSnearestMeas_gain,3);
                
                hDS_interp = real(ifft(XDSnearestMeas_interp));
                % encode to correct rotation. Remove the convert function
                % if the SRIRs are already in N3D normalisation
                hDS_enc = convert_N3D_SN3D(rotateHOA_N3D(convert_N3D_SN3D(hDS_interp,'sn2n'),aziDif,eleDif,0),'n2sn');
                
                % normalise RMS of direct sound
                rmsNormValue = sum(rms(squeeze(nearestMeasDS(:,1,:))) * gains);
                hDS_enc = hDS_enc / rms(hDS_enc(:, 1)) * rmsNormValue;
                hAllDS_interp(1:length_sampDS, :, iInterp) = ...
                    hDS_enc;
            end
            
        case 'minPhase'
            % interpolate direct sound spectrum of the nearest 4 points
            % (if 2d) or 2 points (if 1d)
            nearestMeasDS = ...
                squeeze(hAllDS(1:length_sampDS, 1, idxSorted(1:2^mode_dimension)));
            XDSnearestMeas = fft(nearestMeasDS, nfftDS);
            XDSnearestMeas_smooth = ...
                smoothSpectrum(abs(XDSnearestMeas(1:nfftDS/2+1, :)), fax(:), 3);
            XDSnearestMeas_interp = XDSnearestMeas_smooth * gains;
            hMinPhase = designMinPhase(XDSnearestMeas_interp);
            
            % encode to correct rotation. Remove the convert function if 
            % the SRIRs are in N3D normalisation
            hDS_enc = hMinPhase(1:length_sampDS) * ...
                convert_N3D_SN3D(evalSH(N_interp, [aziTarget, eleTarget]),'n2sn');
            
            % normalise RMS of direct sound
            rmsNormValue = sum(rms(squeeze(nearestMeasDS(:,1,:))) * gains);
            hDS_enc = hDS_enc / rms(hDS_enc(:, 1)) * rmsNormValue;
            hAllDS_interp(1:length_sampDS, :, iInterp) = ...
                hDS_enc;
    end
    
    %% interpolate early reflections
    if sortedDistances(1) == 0 % if it's an exact point, don't interpolate
        hAllER_interp(earlyRefsStart_samp:earlyRefsEnd_samp,:,iInterp) = ...
            earlyRefsIsolated(:,1:numAmbiChannels,1);
    else
        % take fft of nearest measurements, weight by gains
        XERnearestMeas = fft(earlyRefsIsolated);
        XERnearestMeas_gain = XERnearestMeas .* gainsSH;
        
        % then get the weighted nearest measurements in sectors
        XERnearestMeas_gain_secs = ...
            zeros(size(XERnearestMeas_gain,1),numSecs,size(XERnearestMeas_gain,3));
        for i = 1:2^mode_dimension
            XERnearestMeas_gain_secs(:,:,i) = XERnearestMeas_gain(:,:,i) * w_nm;
        end
        
        HER_interp_secs = zeros(size(XERnearestMeas_gain,1),numSecs);
        for secIdx = 1:numSecs % for each sector:
            XERnearestMeas_gain_sec = XERnearestMeas_gain_secs(:, secIdx,:);
            XERnearestMeas_interp_sec = sum(XERnearestMeas_gain_sec,3);
            
            % ERB frequency bands
            ERBfreqs = flip(-(erb_Q*erb_min)+exp((1:n_bandsER)'*(-log(fs/2 + erb_Q*erb_min) + ...
                log(erb_low_freq + erb_Q*erb_min))/n_bandsER) * (fs/2 + erb_Q*erb_min));
            freqs = 1:fs/nfftER:fs;
            freqsInd = zeros(n_bandsER,1);
            for i = 1:length(ERBfreqs)
                freqsInd(i) = find(freqs>=ERBfreqs(i),1);
            end
            freqsInd(1) = 1; freqsInd(end) = find(freqs>=22000,1); % limit the first and last
            
            % over frequency bands, obtain the RMS magnitude difference of the interpolated and nearest IRs
            rmsER_target = zeros(n_bandsER-1,1);
            rmsER_current = zeros(n_bandsER-1,1);
            for j = 1:n_bandsER-1 % calculate rms of freq bands
                rmsER_target(j,:) = sum(rms(XERnearestMeas_gain_sec(freqsInd(j):freqsInd(j+1),1,:)),3);
                rmsER_current(j,:) = rms(XERnearestMeas_interp_sec(freqsInd(j):freqsInd(j+1),1));
            end
            rmsDiffER = rmsER_target ./ rmsER_current;
            
            % create array of RMS equalisation gains
            rmsGainsER = zeros(length(XERnearestMeas_interp_sec),1);
            for j = 1:n_bandsER-1
                if j == n_bandsER-1 % last band
                    rmsGainsER(freqsInd(j):freqsInd(j+1),1) = ...
                        linspace(rmsDiffER(j,1),rmsDiffER(j,1),length(freqsInd(j):freqsInd(j+1)));
                else
                    rmsGainsER(freqsInd(j):freqsInd(j+1),1) = ...
                        linspace(rmsDiffER(j,1),rmsDiffER(j+1,1),length(freqsInd(j):freqsInd(j+1)));
                end
            end
            % fade out the 20khz band to the nfft/2 band with a cosine window
            rmsGainsER(freqsInd(end)+1:nfftER/2) = (cos((0:length(rmsGainsER(freqsInd(end)+1:nfftER/2))-1)' /...
                (length(rmsGainsER(freqsInd(end)+1:nfftER/2))-1) * pi / 2).^2)*rmsDiffER(end);
            rmsGainsER(nfftER/2+1:end) = flip(rmsGainsER(1:nfftER/2));
            
            % apply gains
            HER_interp_sec = XERnearestMeas_interp_sec.*rmsGainsER;
            HER_interp_secs(:,secIdx) = HER_interp_sec;
        end
        
        % back into SH domain (still in Freq domain though)
        HER_interp = HER_interp_secs * w_nm';
        
        % needs some normalisation due to the beamWeightsMaxEV
        w_n_gain = w_n.^2 * length(secDirs(:,1)) ./ ((0:1:N_interp) * 2 + 1)';
        for i = 1:length(HER_interp(1,:))
            n = ceil(sqrt(i)-1);
            HER_interp(:,i) = HER_interp(:,i) / w_n_gain(n+1);
        end
        
        % return to time-domain and window
        hER_interp = real(ifft(HER_interp));
        
        % normalise interpolated early reflections -- for each SH channel
        rmsNormValueER = sum(rms(earlyRefsIsolated) .* repmat(gainsSH,1,numAmbiChannels),3);
        hER_interp = hER_interp ./ rms(hER_interp) .* rmsNormValueER;
        
        hAllER_interp(earlyRefsStart_samp:earlyRefsEnd_samp, :, iInterp) = ...
            hER_interp;
    end
    
    %% Interpolate late reverb
    if sortedDistances(1) == 0 % if it's an exact point, don't interpolate
        hAllLR_interp(lateRevStart_samp:end,:,iInterp) = ...
            lateRevIsolated(:, 1:numAmbiChannels, 1);
    else
        XLRnearestMeas = fft(lateRevIsolated);
        
        % linear interpolation
        XLRnearestMeas_gain = XLRnearestMeas .* gainsSH;
        XLRnearestMeas_interp = sum(XLRnearestMeas_gain,3);
        
        % ERB frequency bands
        ERBfreqs = flip(-(erb_Q*erb_min)+exp((1:n_bandsLR)'*(-log(fs/2 + erb_Q*erb_min) + ...
            log(erb_low_freq + erb_Q*erb_min))/n_bandsLR) * (fs/2 + erb_Q*erb_min));
        freqs = 1:fs/nfftLR:fs;
        freqsInd = zeros(n_bandsLR,1);
        for i = 1:length(ERBfreqs)
            freqsInd(i) = find(freqs>=ERBfreqs(i),1);
        end
        freqsInd(1) = 1; freqsInd(end) = find(freqs>=22000,1); % limit the first and last
        
        % over frequency bands, obtain the RMS magnitude difference of the interpolated and nearest IRs
        rmsLR_target = zeros(n_bandsLR-1,1);
        rmsLR_current = zeros(n_bandsLR-1,1);
        for j = 1:n_bandsLR-1 % calculate rms of freq bands
            rmsLR_target(j,:) = sum(rms(XLRnearestMeas_gain(freqsInd(j):freqsInd(j+1),1,:)),3);
            rmsLR_current(j,:) = rms(XLRnearestMeas_interp(freqsInd(j):freqsInd(j+1),1));
        end
        rmsDiffLR = rmsLR_target ./ rmsLR_current;
        
        % create array of RMS equalisation gains
        rmsGainsLR = zeros(length(XLRnearestMeas_interp),1);
        for j = 1:n_bandsLR-1
            if j == n_bandsLR-1 % last band
                rmsGainsLR(freqsInd(j):freqsInd(j+1),1) = ...
                    linspace(rmsDiffLR(j,1),rmsDiffLR(j,1),length(freqsInd(j):freqsInd(j+1)));
            else
                rmsGainsLR(freqsInd(j):freqsInd(j+1),1) = ...
                    linspace(rmsDiffLR(j,1),rmsDiffLR(j+1,1),length(freqsInd(j):freqsInd(j+1)));
            end
        end
        % fade out the 20khz band to the nfft/2 band with a cosine window
        rmsGainsLR(freqsInd(end)+1:nfftLR/2) = (cos((0:length(rmsGainsLR(freqsInd(end)+1:nfftLR/2))-1)' /...
            (length(rmsGainsLR(freqsInd(end)+1:nfftLR/2))-1) * pi / 2).^2)*rmsDiffLR(end);
        rmsGainsLR(nfftLR/2+1:end) = flip(rmsGainsLR(1:nfftLR/2));
        
        % apply gains, return to time domain and windowing at start and end of LR
        hLR_interp = real(ifft(XLRnearestMeas_interp.*repmat(rmsGainsLR,[1,numAmbiChannels])));
        
        % normalise late reverb for each SH channel
        rmsNormValueLR = sum(rms(lateRevIsolated) .* repmat(gainsSH,1,numAmbiChannels),3);
        hLR_interp = hLR_interp ./ rms(hLR_interp) .* rmsNormValueLR;
        
        hAllLR_interp(lateRevStart_samp:end, :, iInterp) = ...
            hLR_interp;
    end
    
    %% Construct the final IRs
    srirs_interp(:, :, iInterp) = hAllDS_interp(:, :, iInterp) + ...
        hAllER_interp(:, :, iInterp) +...
        hAllLR_interp(:, :, iInterp);
    
end
end