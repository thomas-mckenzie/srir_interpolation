function [srirs_interp,pos_interp] = interpolate_SRIRs_win(srirs_input,pos_input,resolution_new,fs,interp_modeDS, region_limits_db, length_sampDS,fade_sampDS2Regions, fade_sampRegions,N_interp)
% Perceptual interpolation of SRIRs
%   srirs_input - input format [samples, sh_channels, num_measurements]

if nargin<10; N_interp = sqrt(size(srirs_input,2))-1;    end
if nargin<9; fade_sampRegions = fs/100;                   end
if nargin<8; fade_sampDS2Regions = 20;                       end
if nargin<7; length_sampDS = 200;                       end
if nargin<6; region_limits_db = [-10, -20, -30, -40, -50]; end
if nargin<5; interp_modeDS = 'minPhase';                end

%% init

irLength_samp = size(srirs_input,1);
numMeasurements = size(srirs_input,3);
numAmbiChannels = size(srirs_input,2);

nfftDS = 2^nextpow2(length_sampDS);


% ERB parameters (Glasberg and Moore)
erb_Q = 9.26449;
erb_min = 24.7;
erb_low_freq = 10;
n_bands = 48;

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

numRegions = length(region_limits_db);

region_limits_samp = zeros(numRegions, numMeasurements);

[B_decay, A_decay] = octdsgn(1000, fs, 3);
for i = 1:numMeasurements
    % filter and schroeder integration to get EDC. W channel, 1kHz freq band
    srir_input_filt = filter(B_decay,A_decay,srirs_input(:,1,i));
    int_sch = 1/length(srir_input_filt) * cumtrapz(flip(srir_input_filt/max(abs(srir_input_filt))).^2);
    edc = flip(int_sch(2:end));
    edc = edc/max(edc);
    
    for iRegion = 1:numRegions 
    ind_at_t_L = find(db(edc, 'power')<region_limits_db(iRegion),1);
    region_limits_samp(iRegion, i) = ceil(ind_at_t_L/500) * 500; % nearest 500 (empirically chosen)
    end
    
end

region_limits_samp(region_limits_samp>irLength_samp) = irLength_samp;

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

winDS = [ones(length_sampDS-fade_sampDS2Regions-1, 1); cos((0:fade_sampDS2Regions)' / fade_sampDS2Regions * pi / 2).^2];
wins = zeros(irLength_samp,numAmbiChannels,numRegions, numMeasurements);

for i = 1:numMeasurements
    
    for iRegion = 1:numRegions
        if iRegion == 1
            % length of the window
            lenWin = region_limits_samp(1,i) - length_sampDS;
            
            endWin(iRegion, i) = region_limits_samp(1, i) + fade_sampRegions / 2;
            startWin(iRegion, i) = length_sampDS - fade_sampDS2Regions / 2;
            
            % at the boundary, the fade in is at half the length 
            win = [sin((0:fade_sampDS2Regions-1)' / fade_sampDS2Regions * pi /2).^2; ...
            ones(lenWin-fade_sampRegions/2-fade_sampDS2Regions/2, 1); ...
            cos((0:fade_sampRegions-1)' / fade_sampRegions * pi /2).^2];
            

        else 
            lenWin = region_limits_samp(iRegion,i) - region_limits_samp(iRegion-1,i);
            
            endWin(iRegion, i)  = region_limits_samp(iRegion, i) + fade_sampRegions / 2;
            
            startWin(iRegion, i)  = region_limits_samp(iRegion-1, i) - fade_sampRegions / 2;
        
            % at the boundary, the fade in is at half the length 
            win = [sin((0:fade_sampRegions-1)' / fade_sampRegions * pi /2).^2; ...
            ones(lenWin-fade_sampRegions, 1); ...
            cos((0:fade_sampRegions-1)' / fade_sampRegions * pi /2).^2];
        end    
            wins(startWin(iRegion, i):endWin(iRegion, i) -1,:,iRegion, i) = ...
                repmat(win,[1,numAmbiChannels]);  
    end
end

wins(irLength_samp+1:end, :, :, :) = [];
endWin = min(endWin, irLength_samp);

hAllDirectSound = zeros(size(srirs_input, 1), numAmbiChannels, numMeasurements);
hAllDirectSound(1:length(winDS),:,:) = srirs_input(1:length(winDS),1:numAmbiChannels,:) .* repmat(winDS,1,numAmbiChannels,numMeasurements);

%% Spherical Filter Bank / Beamforming
[~, secDirs] = getTdesign(2*N_interp+1);
secDirs = rad2deg(secDirs);
numSecs = size(secDirs, 1);

% Get beam weights. Using max energy vector weights from spherical array processing toolbox
w_n = beamWeightsMaxEV(N_interp);  

w_nm = zeros(size(srirs_input, 2),size(secDirs,1));
for i = 1:length(secDirs(:,1))
    w_nm(:,i) = rotateAxisCoeffs(w_n, pi/2-(deg2rad(secDirs(i,2))), deg2rad(secDirs(i,1)), 'real')';
end

%% Interpolation

fax = 0:fs/nfftDS:fs/2;
hAllDirectSoundInterpolated = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);

srirs_interp = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);

for iInterp = 1:numInterpolatedPoints
    % calculate closest measurements and distances for interpolation
    [sortedDistances, idxSorted] = sort(vecnorm(pos_input - pos_interp(iInterp, :), 2, 2));
    idxNearest = idxSorted(1); % index of closest measurement
    
    gains = 1 ./ (sortedDistances(1:2^mode_dimension)+eps);
    gains = gains / sum(abs(gains));
    gainsSH(1,1,:) = gains; % for the SH channel gains (ER and LR)
    
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
        
    %% interpolate reverb
    
    hAllReverbInterpolated = zeros(irLength_samp, numAmbiChannels, numInterpolatedPoints);
    
for iRegion = 1:numRegions
    
    % closest points
    idxSelected =  idxSorted(1:2^mode_dimension);
    
    selectedStartWin = min(startWin(iRegion, idxSelected));
    selectedEndWin = max(endWin(iRegion, idxSelected));
    
    regionLength = selectedEndWin-selectedStartWin+1;
   
    % power of 2 for computational speed, one more power to account for some time aliasing
    nfft = 2^(nextpow2(regionLength)+1); 
    
    hRegion = srirs_input(selectedStartWin:selectedEndWin, ...
        1:numAmbiChannels, idxSelected).* ...
        squeeze(wins(selectedStartWin:selectedEndWin, :, iRegion, idxSelected));

    if sortedDistances(1) == 0 % if it's an exact point, don't bother interpolating
        hAllReverbInterpolated(selectedStartWin:selectedEndWin, : ,iInterp) = ...
            hAllReverbInterpolated(selectedStartWin:selectedEndWin,:,iInterp)+ ...
            hRegion(:, :, 1);
    else
        
        % take fft of nearest measurements, weight by gains
        XRefsNearestPoints = fft(hRegion, nfft);
        XRefsNearestPoints_gain = XRefsNearestPoints .* gainsSH;
        
        % then get the weighted nearest measurements in sectors
        XRefsNearestPoints_gain_secs = zeros(length(XRefsNearestPoints_gain(:,1,1)),numSecs,length(XRefsNearestPoints_gain(1,1,:)));
        for i = 1:2^mode_dimension
            XRefsNearestPoints_gain_secs(:,:,i) = XRefsNearestPoints_gain(:,:,i) * w_nm;
        end
        
        Hinterp_secs = zeros(length(XRefsNearestPoints_gain(:,1,1)),numSecs);
        
        for secIdx = 1:numSecs % for each sector
            
            XnearestMeas_gain_sec = XRefsNearestPoints_gain_secs(:, secIdx,:);
            XnearestMeas_interp_sec = sum(XnearestMeas_gain_sec,3);
            
            % ERB frequency bands
            ERBfreqs = flip(-(erb_Q*erb_min)+exp((1:n_bands)'*(-log(fs/2 + erb_Q*erb_min) + ...
                log(erb_low_freq + erb_Q*erb_min))/n_bands) * (fs/2 + erb_Q*erb_min));
            freqs = 1:fs/nfft:fs;
            freqsInd = zeros(n_bands,1);
            for i = 1:length(ERBfreqs)
                freqsInd(i) = find(freqs>=ERBfreqs(i),1);
            end
            freqsInd(1) = 1; freqsInd(end) = find(freqs>=22000,1); % limit the first and last
            
            % over frequency bands, obtain the RMS magnitude difference of the interpolated and nearest IRs
            rmstarget = zeros(n_bands-1,1);
            rmscurrent = zeros(n_bands-1,1);
            for j = 1:n_bands-1 % calculate rms of freq bands
                rmstarget(j,:) = sum(rms(XnearestMeas_gain_sec(freqsInd(j):freqsInd(j+1),1,:)),3);
                rmscurrent(j,:) = rms(XnearestMeas_interp_sec(freqsInd(j):freqsInd(j+1),1));
            end
            rmsDiff = rmstarget ./ rmscurrent;
            
            % create array of RMS equalisation gains
            rmsGains = zeros(length(XnearestMeas_interp_sec),1);
            for j = 1:n_bands-1
                if j == n_bands-1 % last band
                    rmsGains(freqsInd(j):freqsInd(j+1),1) = ...
                        linspace(rmsDiff(j,1),rmsDiff(j,1),length(freqsInd(j):freqsInd(j+1)));
                else
                    rmsGains(freqsInd(j):freqsInd(j+1),1) = ...
                        linspace(rmsDiff(j,1),rmsDiff(j+1,1),length(freqsInd(j):freqsInd(j+1)));
                end
            end
            % fade out the 20khz band to the nfft/2 band with a cosine window
            rmsGains(freqsInd(end)+1:nfft/2) = (cos((0:length(rmsGains(freqsInd(end)+1:nfft/2))-1)' /...
                (length(rmsGains(freqsInd(end)+1:nfft/2))-1) * pi / 2).^2)*rmsDiff(end);
            rmsGains(nfft/2+1:end) = flip(rmsGains(1:nfft/2));
            
            % apply gains
            Hinterp_sec = XnearestMeas_interp_sec.* rmsGains;
            Hinterp_secs(:,secIdx) = Hinterp_sec;
        end
        
        % back into SH domain (still in Freq domain though)
        Hinterp = Hinterp_secs * w_nm';
        
        % needs some normalisation due to the beamWeightsMaxEV
        w_n_gain = w_n.^2 * length(secDirs(:,1)) ./ ((0:1:N_interp) * 2 + 1)';
        for i = 1:length(Hinterp(1,:))
            n = ceil(sqrt(i)-1);
            Hinterp(:,i) = Hinterp(:,i) / w_n_gain(n+1);
        end
        
        % return to time-domain and window
        hInterp = real(ifft(Hinterp));
        
        % normalise interpolated early reflections -- for each SH channel
        rmsNormValue = sum(rms(hRegion) .* repmat(gainsSH,1,numAmbiChannels),3);
        hInterp = hInterp ./ rms(hInterp) .* rmsNormValue;
        
        hAllReverbInterpolated(selectedStartWin:selectedEndWin, : ,iInterp) = ...
            hAllReverbInterpolated(selectedStartWin:selectedEndWin,:,iInterp) + hInterp(1:regionLength, :);
        
        figure(10), plot(hAllReverbInterpolated(:, 1, iInterp))
        pause(.1)
        iInterp
    end    
end
    %% Construct the final IRs
    % add three parts together
    srirs_interp(:, :, iInterp) = hAllDirectSoundInterpolated(:, :, iInterp) + ...
        hAllReverbInterpolated(:, :, iInterp);
end