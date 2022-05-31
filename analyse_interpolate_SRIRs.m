% analyse interpolated IRs

    directSoundLength_samp = 200;
    earlyRefsCutoff_samp = 10000;
    N_interp = 4;

figure;
title(['RMS and DRR values. Measurement reduction = ',num2str(measReduction),'.']);


subplot(2,4,1)
plot(rms(abs(squeeze(hAllInterpolated(1:directSoundLength_samp,1,:)))));
hold on
plot(rms(abs(squeeze(hCorrect(1:directSoundLength_samp,1,:)))));
title('Direct sound RMS W channel')

subplot(2,4,2)
plot(rms(abs(squeeze(hAllInterpolated(directSoundLength_samp:earlyRefsCutoff_samp,1,:)))));
hold on
plot(rms(abs(squeeze(hCorrect(directSoundLength_samp:earlyRefsCutoff_samp,1,:)))));
title('Early Reflections RMS W channel')

erEnd = earlyRefsCutoff_samp+directSoundLength_samp;
subplot(2,4,3)
plot(rms((squeeze(hAllInterpolated(erEnd:end,1,:)))));
hold on
plot(rms((squeeze(hCorrect(erEnd:end,1,:)))));
title('Late Reverb RMS W channel')

subplot(2,4,4)
plot(rms(abs(squeeze(hAllInterpolated(1:directSoundLength_samp,1,:))))./rms(abs(squeeze(hAllInterpolated(directSoundLength_samp:end,1,:)))));
hold on;plot(rms(abs(squeeze(hCorrect(1:directSoundLength_samp,1,:))))./rms(abs(squeeze(hCorrect(directSoundLength_samp:end,1,:)))));
title('DRR W channel')


subplot(2,4,5)
plot(squeeze(mean(rms(hAllInterpolated(1:directSoundLength_samp,:,:)),2)));
hold on
plot(squeeze(mean(rms(hCorrect(1:directSoundLength_samp,:,:)),2)));
title('Direct sound RMS All channels')

subplot(2,4,6)
plot(squeeze(mean(rms(hAllInterpolated(directSoundLength_samp:earlyRefsCutoff_samp,:,:)),2)));
hold on
plot(squeeze(mean(rms(hCorrect(directSoundLength_samp:earlyRefsCutoff_samp,:,:)),2)));
title('Early Reflections RMS All channels')

subplot(2,4,7)
plot(squeeze(mean(rms(hAllInterpolated(erEnd:end,:,:)),2)));
hold on
plot(squeeze(mean(rms(hCorrect(erEnd:end,:,:)),2)));
title('Late Reverb RMS All channels')


subplot(2,4,8)
plot(mean(squeeze(rms(abs(squeeze(hAllInterpolated(1:directSoundLength_samp,:,:))))./rms(abs(squeeze(hAllInterpolated(directSoundLength_samp:end,:,:)))))));
hold on;plot(mean(squeeze(rms(abs(squeeze(hCorrect(1:directSoundLength_samp,:,:))))./rms(abs(squeeze(hCorrect(directSoundLength_samp:end,:,:)))))));
title('DRR All channels')

%% time domain plots
% channelToPlot = lastChannel; % good to check not just W channel (1) to see if the higher order channels are good too
channelToPlot = 1; % good to check not just W channel (1) to see if the higher order channels are good too
plot_length = 300; % plot length in ms

plot_length_samples = plot_length / 1000 * fs;
[X,Y] = meshgrid(-250:5:250,1/48:1/48:plot_length);
irTrunc = hCorrect;
figure;subplot(1,2,1)
surf(X,Y,real(10*log10(abs(squeeze(irTrunc(1:plot_length_samples,channelToPlot,:))))),'EdgeColor','none');

ylabel('Time (ms)');
c = flip(bone);caxis([-40 0]);colormap(c);h = colorbar;ylabel(h, 'Amplitude (dB)');
view(0,90);
set(gcf, 'Color', 'w');pbaspect([1.7 1 1]);set(gca, 'fontsize', 12);
title(['time-domain - original']);

% interpolated
irTrunc = hAllInterpolated;
subplot(1,2,2)
surf(X,Y,real(10*log10(abs(squeeze(irTrunc(1:plot_length_samples,channelToPlot,:))))),'EdgeColor','none');

ylabel('Time (ms)');
c = flip(bone);caxis([-40 0]);colormap(c);h = colorbar;ylabel(h, 'Amplitude (dB)');
view(0,90);
set(gcf, 'Color', 'w');pbaspect([1.7 1 1]);set(gca, 'fontsize', 12);
title(['time-domain - interpolated. Measurement reduction = ',num2str(measReduction),'.']);

%% DOA

% original
irTrunc = hCorrect(:,1:(N_interp+1)^2,:);

% Tuneable parameters:
degreeResolution = 5;
order = N_interp;
nSrc = 7;
numSamps = 4800*3; % 0.3 seconds
highPassFilterFreq = 3000;
kappa = 40;

grid_dirs = grid2dirs(degreeResolution,degreeResolution,0,0); % Grid of directions to evaluate DoA estimation
P_src = diag(ones(numSamps,1));
[~,filtHi,~] = ambisonic_crossover(highPassFilterFreq,fs);

doa_est = zeros(nSrc,2,length(irTrunc(1,1,:)));
doa_est_P = zeros(nSrc,length(irTrunc(1,1,:)));

for i = 1:length(irTrunc(1,1,:))
    Y_src = filter(filtHi,1,irTrunc(1:numSamps,:,i)); % high pass filter
    
    stVec = Y_src';
    sphCOV = stVec*P_src*stVec' + 1*eye((order+1)^2)/(4*pi);
    
    % DoA estimation
    [~, est_dirs_pwd,est_dirs_P] = sphPWDmap(sphCOV, grid_dirs, nSrc,kappa);
    
    % convert to degs from rads
    est_dirs_pwd = est_dirs_pwd*180/pi;
    
    % flip -ve values near -180 to +ve values
    negativeFlipLimit = -170;
    for j = 1:length(est_dirs_pwd(:,1))
        for k = 1:length(est_dirs_pwd(1,:))
            if est_dirs_pwd(j,k) < negativeFlipLimit
                est_dirs_pwd(j,k) = est_dirs_pwd(j,k) + 360;
            end
        end
    end
    doa_est(:,:,i) = est_dirs_pwd;
    doa_est_P(:,i) = est_dirs_P;
end

normalized_doa_est_P = doa_est_P ./ max(doa_est_P,[],1);
normalized_doa_est_P_dB = mag2db(normalized_doa_est_P)/2;

plot_thresh = -10;
normalized_doa_est_P_dB( normalized_doa_est_P_dB < plot_thresh ) = plot_thresh;
c_truncation = 20;
cmap = flip(parula(256));
c = cmap(c_truncation:end,:); % truncate yellow

doa_01 = rescale(normalized_doa_est_P_dB, 'InputMin',plot_thresh);
cspace = linspace(0,1,size(c,1));
doa_color(:,:,1) = interp1(cspace, c(:,1), doa_01);
doa_color(:,:,2) = interp1(cspace, c(:,2), doa_01);
doa_color(:,:,3) = interp1(cspace, c(:,3), doa_01);

   h= figure;subplot(2,1,1)
for i = nSrc:-1:1
    s = scatter(-250:5:250,squeeze(doa_est(i,1,:)),25,...
        squeeze(doa_color(i,:,:)),'filled');
    hold on
end

ylabel('Azimuth (°)');
xlabel('<-- Storage        Measurement position (cm)        Stairwell -->');

ylim([negativeFlipLimit (negativeFlipLimit+360)]);yticks(-180:45:180);
xlim([-250 250]);xticks(-250:50:250);colormap(c);
k = colorbar;ylabel(k,'Normalised power (dB)');caxis([plot_thresh 0]);
set(gcf, 'Color', 'w');pbaspect([1.7 1 1]);
box on;grid on;set(gca,'FontSize',16)
title('DoA - original');

% Interpolated
irTrunc = hAllInterpolated;
doa_est = zeros(nSrc,2,length(irTrunc(1,1,:)));
doa_est_P = zeros(nSrc,length(irTrunc(1,1,:)));
for i = 1:length(irTrunc(1,1,:))
    Y_src = filter(filtHi,1,irTrunc(1:numSamps,:,i)); % high pass filter
    
    stVec = Y_src';
    sphCOV = stVec*P_src*stVec' + 1*eye((order+1)^2)/(4*pi);
    
    % DoA estimation
    [~, est_dirs_pwd,est_dirs_P] = sphPWDmap(sphCOV, grid_dirs, nSrc,kappa);
    
    % convert to degs from rads
    est_dirs_pwd = est_dirs_pwd*180/pi;
    
    % flip -ve values near -180 to +ve values
    negativeFlipLimit = -170;
    for j = 1:length(est_dirs_pwd(:,1))
        for k = 1:length(est_dirs_pwd(1,:))
            if est_dirs_pwd(j,k) < negativeFlipLimit
                est_dirs_pwd(j,k) = est_dirs_pwd(j,k) + 360;
            end
        end
    end
    doa_est(:,:,i) = est_dirs_pwd;
    doa_est_P(:,i) = est_dirs_P;
end

normalized_doa_est_P = doa_est_P ./ max(doa_est_P,[],1);
normalized_doa_est_P_dB = mag2db(normalized_doa_est_P)/2;

normalized_doa_est_P_dB( normalized_doa_est_P_dB < plot_thresh ) = plot_thresh;

doa_01 = rescale(normalized_doa_est_P_dB, 'InputMin',plot_thresh);

cspace = linspace(0,1,size(c,1));
doa_color(:,:,1) = interp1(cspace, c(:,1), doa_01);
doa_color(:,:,2) = interp1(cspace, c(:,2), doa_01);
doa_color(:,:,3) = interp1(cspace, c(:,3), doa_01);

subplot(2,1,2)
for i = nSrc:-1:1
    s = scatter(-250:5:250,squeeze(doa_est(i,1,:)),25,...
        squeeze(doa_color(i,:,:)),'filled');
    hold on
end

ylabel('Azimuth (°)');
xlabel('<-- Storage        Measurement position (cm)        Stairwell -->');
ylim([negativeFlipLimit (negativeFlipLimit+360)]);yticks(-180:45:180);
xlim([-250 250]);xticks(-250:50:250);colormap(c);
k = colorbar;ylabel(k,'Normalised power (dB)');caxis([plot_thresh 0]);
set(gcf, 'Color', 'w');pbaspect([1.7 1 1]);
box on;grid on;set(gca,'FontSize',16)

title(['DoA - interpolated. Measurement reduction = ',num2str(measReduction),'.']);
