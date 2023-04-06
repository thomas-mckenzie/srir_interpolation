function [srirs_interp,pos_interp_src,dir_interp_src] = interpolate_SRIRs_source(srirs_input,pos_input,division,fs,interp_modeDS,length_sampDS,fade_sampDS2ER,fade_sampER2LR,N_interp)
% Perceptually informed interpolation of SRIRs with different source
% positions
% 
% interpolate_SRIRs_source()
% srirs_input - input format [samples, sh_channels, num_measurements]
% 
% The user is directed to the following paper: 
% McKenzie, T. and Schlecht, S. J. "Source position interpolation of 
%     spatial room impulse responses" AES 154th Convention, Espoo, Helsinki, 
%     Finland, 2023. 
% 
% Requires Spherical Array Processing Toolbox (Politis, 2015)
% (https://github.com/polarch/Spherical-Array-Processing)
% 
% Thomas McKenzie, University of Edinburgh, 2023. thomas.mckenzie@ed.ac.uk

if nargin<9; N_interp = sqrt(size(srirs_input,2))-1;    end
if nargin<8; fade_sampER2LR = fs/100;                   end
if nargin<7; fade_sampDS2ER = 20;                       end
if nargin<6; length_sampDS = 200;                       end
if nargin<5; interp_modeDS = 'meanSpectrum';            end

reinsertTimeAlignmentDelayAfterInterp = 0;

%% init
irLength_samp = size(srirs_input,1);
numMeas = size(srirs_input,4); % 4th dimension for num of interpolation measurements
numRecPosns = size(srirs_input,3); % 3rd dimension for number of receiver positions
numAmbiChans = (N_interp+1)^2; % number of Ambisonic channels
nfftDS = 2^nextpow2(length_sampDS);

% ERB parameters (Glasberg and Moore)
erb_Q = 9.26449;erb_min = 24.7;erb_low_freq = 10;
n_bandsDS = 48;n_bandsER = 48;n_bandsLR = 48;

%% Interpolated positions
mode_dimension = 1; % fixed for now

% for the division: interval is how many RIRs are we interpolating here
interval = (pos_input(1,:) - pos_input(2,:))/(division);
for i = 1:3
    if interval(i) == 0
        pos_interp_src(:,i) = pos_input(1,i);
    else
        pos_interp_src(:,i) = pos_input(1,i) : - interval(i) : pos_input(2,i);
    end
end
numInterpPts = size(pos_interp_src,1);

%% Reorder input SRIRs
srirs_input = permute(srirs_input,[1 2 4 3]); % reorder here: indexes are now [samples, SH channels, source positions, receiver positions]
srirs_interp = zeros(irLength_samp, numAmbiChans, numInterpPts,numRecPosns);

%% Time align onsets
initialInputGap = 50; % pre-delay before the detected onset for IRs
onsetThreshold = 0.1; % empirically chosen threshold used for truncation of IRs, chosen by visual observation of onset and noise floor amplitudes

indexMin = zeros(size(srirs_input,3),size(srirs_input,4));
for i = 1:size(srirs_input,3)
    for j = 1:size(srirs_input,4)
        indexMin(i,j) = find(abs(srirs_input(:,1,i,j)/max(abs(srirs_input(:,1,i,j))))>onsetThreshold,1);
        if indexMin(i,j) <= initialInputGap
            initialInputGap = 0; % catch it if the onset is already below the input gap.
        end
        srirs_input(:,:,i,j) = [srirs_input(indexMin(i,j) - initialInputGap:end,:,i,j) ;  zeros((indexMin(i,j) - initialInputGap-1),numAmbiChans)];
    end
end

%% Interpolate between source positions for each receiver position
for iRec = 1:numRecPosns

    %% Calculate the cutoff between early refs and late reverb
    % based on when the energy decay curve amplitude is lower than a target

    ER_target = 0.1; % empirically chosen value. could be 1 / exp(1);
    ERCutoff_samp = zeros(1,numMeas);
    [B_decay, A_decay] = octave_filter(1000, fs, 3);
    for i = 1:numMeas
        % filter and schroeder integration to get EDC. W channel, 1kHz freq band
        srir_input_filt = filter(B_decay,A_decay,srirs_input(:,1,i,iRec));
        int_sch = 1/length(srir_input_filt) * cumtrapz(flip(srir_input_filt/max(abs(srir_input_filt))).^2);
        edc = flip(int_sch(2:end));
        edc = edc/max(edc);

        ind_at_t_L = find(edc<ER_target,1);
        ERCutoff_samp(i) = ceil(ind_at_t_L/1000) * 1000; % nearest 1000 (empirically chosen)
    end

    %% input SRIR DoA (for interpolating direct sound)

    directSoundDirections = zeros(numMeas, 3);
    for iMeas = 1:numMeas
        srir_doa = srirs_input(1:length_sampDS, 1:4, iMeas, iRec);

        % intensity vector DoA analysis
        doa_samp = srir_doa(:, 1) .* [srir_doa(:, 4), srir_doa(:, 2), srir_doa(:, 3)];
        idxNonZero = find(sum(doa_samp ~= 0, 2));
        doa = mean(doa_samp(idxNonZero, :));
        doa = doa ./ vecnorm(doa, 2, 2);
        directSoundDirections(iMeas, :) = doa;
    end

    %% cut out direct sound and make windows for early refs and late reverberation
    % windows for ER and LR can vary with the ER window size (earlyRefsCutoff_samp)

    winDS = [ones(length_sampDS-fade_sampDS2ER-1, 1); cos((0:fade_sampDS2ER)' / fade_sampDS2ER * pi / 2).^2];
    winER = zeros(max(ERCutoff_samp)+fade_sampER2LR+length_sampDS-fade_sampDS2ER,numAmbiChans,numMeas);
    winLR = ones(irLength_samp,numAmbiChans,numMeas);

    winERlength = ERCutoff_samp+fade_sampER2LR+length_sampDS-fade_sampDS2ER;
    for i = 1:numMeas
        winER(1:winERlength(i),:,i) = repmat([zeros(length_sampDS-fade_sampDS2ER-1,1);...
            sin((0:fade_sampDS2ER)' / fade_sampDS2ER * pi / 2).^2; ...
            ones(ERCutoff_samp(i)-fade_sampDS2ER-1, 1 ); ...
            cos((0:fade_sampER2LR)' / fade_sampER2LR * pi / 2).^2],[1,numAmbiChans]);

        winLR(1:winERlength(i),:,i) = repmat([zeros(length_sampDS-fade_sampDS2ER-1,1);...
            zeros(ERCutoff_samp(i),1);...
            (1-cos((0:fade_sampER2LR)' / fade_sampER2LR * pi / 2).^2)],[1,numAmbiChans]);
    end

    hAllDS = zeros(size(srirs_input, 1), numAmbiChans, numMeas);
    hAllDS(1:length(winDS),:,:) = squeeze(srirs_input(1:length(winDS),1:numAmbiChans,:,iRec)) .*...
        repmat(winDS,1,numAmbiChans,numMeas);

    %% Spherical Filter Bank / Beamforming
    [~, secDirs] = getTdesign(2*N_interp+1); % needs spherical array processing toolbox
    secDirs = rad2deg(secDirs);
    numSecs = size(secDirs, 1);

    % Get beam weights. Here using max energy vector weights, but can change
    w_n = beamWeightsMaxEV(N_interp); % needs spherical array processing toolbox

    w_nm = zeros(numAmbiChans,numSecs);
    for i = 1:length(secDirs(:,1))
        w_nm(:,i) = rotateAxisCoeffs(w_n, pi/2-(deg2rad(secDirs(i,2))), deg2rad(secDirs(i,1)), 'real')';
    end

    %% Interpolation
    fax = 0:fs/nfftDS:fs/2;
    hAllDS_interp = zeros(irLength_samp, numAmbiChans, numInterpPts);
    hAllER_interp = zeros(irLength_samp, numAmbiChans, numInterpPts);
    hAllLR_interp = zeros(irLength_samp, numAmbiChans, numInterpPts);
    srirs_interp_1rec = zeros(irLength_samp, numAmbiChans, numInterpPts);

    for iInterp = 1:numInterpPts
        % calculate closest measurements and distances for interpolation
        [sortedDistances, idxSorted] = sort(vecnorm(pos_input - pos_interp_src(iInterp, :), 2, 2));
        idxNearest = idxSorted(1); % index of closest measurement

        gains = 1 ./ (sortedDistances(1:2^mode_dimension)+eps);
        gains = gains / sum(abs(gains));
        gainsSH(1,1,:) = gains; % for the SH channel gains (ER and LR)

        % get start and end samples and lengths of windows for RIR portioning
        earlyRefsStart_samp = length_sampDS-fade_sampDS2ER+1;
        earlyRefsEnd_samp = ERCutoff_samp(idxNearest)+fade_sampER2LR+earlyRefsStart_samp-1;
        lateRevStart_samp = ERCutoff_samp(idxNearest)+earlyRefsStart_samp;

        % fft window sizes
        nfftER = length(earlyRefsStart_samp:earlyRefsEnd_samp);
        nfftLR = length(lateRevStart_samp:irLength_samp);

        % isolate early and late reverb, and window
        earlyRefsIsolated = squeeze(srirs_input(earlyRefsStart_samp:earlyRefsEnd_samp, 1:numAmbiChans, idxSorted(1:2^mode_dimension), iRec)).* ...
            repmat(winER(earlyRefsStart_samp:earlyRefsEnd_samp,1:numAmbiChans,idxNearest),[1,1,2^mode_dimension]);
        lateRevIsolated = squeeze(srirs_input(lateRevStart_samp:end, 1:numAmbiChans, idxSorted(1:2^mode_dimension), iRec)).* ...
            repmat(winLR(lateRevStart_samp:end,1:numAmbiChans,idxNearest),[1,1,2^mode_dimension]);

        %% interpolate direct sound
        % DoA interp for direct sound interpolation
        dsDir = zeros(2^mode_dimension,2); % direction of arrivals of the direct sounds
        dsDiff = zeros(2^mode_dimension,2); % differences between the direct sounds and weighted average directional target
        for i = 1:2^mode_dimension
            [dsDir(i,1),dsDir(i,2),~] = cart2sph(directSoundDirections(idxSorted(i),1), ...
                directSoundDirections(idxSorted(i),2), ...
                directSoundDirections(idxSorted(i),3));

            % 2*pi to wrap round circle for always +ve numbers
            if dsDir(i,1) < 0; dsDir(i,1) = dsDir(i,1) + 2*pi; end
        end
        dsTarget = gains' * dsDir;
        dir_interp_src(:,iInterp) = dsTarget;
        for i = 1:2^mode_dimension
            dsDiff(i,:) = angdiff(dsTarget, dsDir(i,:));
        end

        switch interp_modeDS
            case 'meanSpectrum'
                % avg of FFTs
                if sortedDistances(1) == 0 % if it's the point of a measurement, bypass interpolation
                    hAllDS_interp(1:length_sampDS, :, iInterp) = ...
                        hAllDS(1:length_sampDS, 1:numAmbiChans, idxNearest);
                else
                    nearestMeasDS = hAllDS(1:length_sampDS, :, idxSorted(1:2^mode_dimension));

                    % Rotate input SRIRs to interpolated DoA
                    for i = 1:2^mode_dimension % rotate each SRIR separately
                        nearestMeasDS(:,:,i) = rotateHOA_N3D(nearestMeasDS(:,:,i),...
                            rad2deg(dsDiff(i,1)),rad2deg(dsDiff(i,2)),0);
                    end

                    XDSnearestMeas = fft(nearestMeasDS,nfftDS);
                    XDSnearestMeas_gain = XDSnearestMeas .* gainsSH;
                    XDSnearestMeas_interp = sum(XDSnearestMeas_gain,3);

                    % ERB frequency bands
                    ERBfreqs = flip(-(erb_Q*erb_min)+exp((1:n_bandsDS)'*(-log(fs/2 + erb_Q*erb_min) + ...
                        log(erb_low_freq + erb_Q*erb_min))/n_bandsDS) * (fs/2 + erb_Q*erb_min));
                    freqs = 1:fs/nfftDS:fs;
                    freqsInd = zeros(n_bandsDS,1);
                    for i = 1:length(ERBfreqs)
                        freqsInd(i) = find(freqs>=ERBfreqs(i),1);
                    end
                    freqsInd(1) = 1; freqsInd(end) = find(freqs>=22000,1); % limit the first and last

                    % over frequency bands, obtain the RMS magnitude difference of the interpolated and nearest IRs
                    rmsDS_target = zeros(n_bandsDS-1,1);
                    rmsDS_current = zeros(n_bandsDS-1,1);
                    for j = 1:n_bandsDS-1 % calculate rms of freq bands
                        rmsDS_target(j,:) = sum(rms(XDSnearestMeas_gain(freqsInd(j):freqsInd(j+1),1,:)),3);
                        rmsDS_current(j,:) = rms(XDSnearestMeas_interp(freqsInd(j):freqsInd(j+1),1));
                    end
                    rmsDiffDS = rmsDS_target ./ rmsDS_current;

                    % create array of RMS equalisation gains
                    rmsGainsDS = zeros(length(XDSnearestMeas_interp),1);
                    for j = 1:n_bandsDS-1
                        if j == n_bandsDS-1 % last band
                            rmsGainsDS(freqsInd(j):freqsInd(j+1),1) = ...
                                repmat(rmsDiffDS(j,1),1,length(freqsInd(j):freqsInd(j+1)));
                        else
                            rmsGainsDS(freqsInd(j):freqsInd(j+1),1) = ...
                                linspace(rmsDiffDS(j,1),rmsDiffDS(j+1,1),length(freqsInd(j):freqsInd(j+1)));
                        end
                    end
                    % continue the 20khz band to the nfft/2 band and flip 2nd half
                    rmsGainsDS(freqsInd(end)+1:nfftDS/2) = rmsDiffDS(end);
                    rmsGainsDS(nfftDS/2+1:end) = flip(rmsGainsDS(1:nfftDS/2));

                    % apply gains, return to time domain and windowing at start and end of DS
                    hDS_enc = real(ifft(XDSnearestMeas_interp.*repmat(rmsGainsDS,[1,numAmbiChans])));
                    hDS_enc = hDS_enc(1:length_sampDS,:);

                    % normalise RMS of direct sound
                    rmsNormValue = sum(rms(squeeze(nearestMeasDS(:,1,:))) * gains);
                    hDS_enc = hDS_enc / rms(hDS_enc(:, 1)) * rmsNormValue;
                    hAllDS_interp(1:length_sampDS, :, iInterp) = hDS_enc;
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
                    convert_N3D_SN3D(evalSH(N_interp, dsTarget),'n2sn');

                % normalise RMS of direct sound
                rmsNormValue = sum(rms(squeeze(nearestMeasDS(:,1,:))) * gains);
                hDS_enc = hDS_enc / rms(hDS_enc(:, 1)) * rmsNormValue;
                hAllDS_interp(1:length_sampDS, :, iInterp) = hDS_enc;
        end

        %% interpolate early reflections
        if sortedDistances(1) == 0 % if it's an exact point, don't interpolate
            hAllER_interp(earlyRefsStart_samp:earlyRefsEnd_samp,:,iInterp) = ...
                earlyRefsIsolated(:,1:numAmbiChans,1);
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
                            repmat(rmsDiffER(j,1),1,length(freqsInd(j):freqsInd(j+1)));
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
            rmsNormValueER = sum(rms(earlyRefsIsolated) .* repmat(gainsSH,1,numAmbiChans),3);
            hER_interp = hER_interp ./ rms(hER_interp) .* rmsNormValueER;

            hAllER_interp(earlyRefsStart_samp:earlyRefsEnd_samp, :, iInterp) = ...
                hER_interp;
        end

        %% Interpolate late reverb
        if sortedDistances(1) == 0 % if it's an exact point, don't interpolate
            hAllLR_interp(lateRevStart_samp:end,:,iInterp) = ...
                lateRevIsolated(:, 1:numAmbiChans, 1);
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
                        repmat(rmsDiffLR(j,1),1,length(freqsInd(j):freqsInd(j+1)));
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
            hLR_interp = real(ifft(XLRnearestMeas_interp.*repmat(rmsGainsLR,[1,numAmbiChans])));

            % normalise late reverb for each SH channel
            rmsNormValueLR = sum(rms(lateRevIsolated) .* repmat(gainsSH,1,numAmbiChans),3);
            hLR_interp = hLR_interp ./ rms(hLR_interp) .* rmsNormValueLR;

            hAllLR_interp(lateRevStart_samp:end, :, iInterp) = hLR_interp;
        end

        %% Construct the final IRs
        srirs_interp_1rec(:, :, iInterp) = hAllDS_interp(:, :, iInterp) + ...
            hAllER_interp(:, :, iInterp) + hAllLR_interp(:, :, iInterp);

    end
    srirs_interp(:,:,:,iRec) = srirs_interp_1rec;
end

%% reinsert delay from time-alignment?
if reinsertTimeAlignmentDelayAfterInterp
    for j = 1:size(srirs_interp,4)
        delay_interp = round(linspace(indexMin(1,j),indexMin(2,j),size(srirs_interp,3)));
        for i = 1:size(srirs_interp,3)
            srirs_interp(:,:,i,j) = [zeros(delay_interp(i),numAmbiChans) ; srirs_interp(1: end-delay_interp(i),:,i,j)];
        end
    end
end

%% Change ordering back
srirs_interp = permute(srirs_interp,[1 2 4 3]);

end

%% Extra functions
function [B,A] = octave_filter(fc,fs,N)
% Requires the Signal Processing Toolbox. 
% Based on the filter design by Christophe Couvreur, Faculte Polytechnique de Mons (Belgium)
% Designs Butterworth 2Nth-order octave filter based on a bilinear transformation (ANSI S1.1-1986)

if nargin < 3; N = 3; end
if (fc > 0.7*(fs/2)); error('Design not possible'); end

b = pi/2/N/sin(pi/2/N); 
a = (1+sqrt(1+8*b^2))/4/b;
W1 = fc/(fs/2)*sqrt(1/2)/a; 
W2 = fc/(fs/2)*sqrt(2)*a;
[B,A] = butter(N,[W1,W2]); 
end


