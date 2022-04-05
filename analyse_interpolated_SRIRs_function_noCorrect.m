% analyse interpolated IRs
% %{

% earlyRefsCutoff_samp = zeros(1,numMeasurements);
% for i = 1:numMeasurements
%     decay = rir2decay(srirs_input(:,1,i), fs, 1000, 1, 1, 1);
%     ind_at_t_L = find(decay<ER_target,1);
    directSoundLength_samp = 200;
    earlyRefsCutoff_samp = 10000;
    


figure;
title(['RMS and DRR values.'])% Measurement reduction = ',num2str(measReduction),'. interpType = ',interpType]);

t_interp = pos_interp(:,1);
t_input = pos_input(:,1);

subplot(2,4,1)
% figure;
plot(t_interp,rms(abs(squeeze(srirs_interp(1:directSoundLength_samp,1,:)))));
hold on
plot(t_input,rms(abs(squeeze(srirs_input(1:directSoundLength_samp,1,:)))));
title('Direct sound RMS W channel')
%}
subplot(2,4,2)
% figure;
plot(t_interp,rms(abs(squeeze(srirs_interp(directSoundLength_samp:earlyRefsCutoff_samp,1,:)))));
hold on
plot(t_input,rms(abs(squeeze(srirs_input(directSoundLength_samp:earlyRefsCutoff_samp,1,:)))));
title('Early Reflections RMS W channel')

% %{
erEnd = earlyRefsCutoff_samp+directSoundLength_samp;
subplot(2,4,3)
% figure;
plot(t_interp,rms((squeeze(srirs_interp(erEnd:end,1,:)))));
hold on
plot(t_input,rms((squeeze(srirs_input(erEnd:end,1,:)))));
title('Late Reverb RMS W channel')

subplot(2,4,4)
% figure;
plot(t_interp,rms(abs(squeeze(srirs_interp(1:directSoundLength_samp,1,:))))./rms(abs(squeeze(srirs_interp(directSoundLength_samp:end,1,:)))));
hold on;plot(t_input,rms(abs(squeeze(srirs_input(1:directSoundLength_samp,1,:))))./rms(abs(squeeze(srirs_input(directSoundLength_samp:end,1,:)))));
title('DRR W channel')


% %{
subplot(2,4,5)

% figure;
plot(t_interp,squeeze(mean(rms(srirs_interp(1:directSoundLength_samp,:,:)),2)));
hold on
plot(t_input,squeeze(mean(rms(srirs_input(1:directSoundLength_samp,:,:)),2)));
title('Direct sound RMS All channels')
%}
subplot(2,4,6)
% figure;
plot(t_interp,squeeze(mean(rms(srirs_interp(directSoundLength_samp:earlyRefsCutoff_samp,:,:)),2)));
hold on
plot(t_input,squeeze(mean(rms(srirs_input(directSoundLength_samp:earlyRefsCutoff_samp,:,:)),2)));
title('Early Reflections RMS All channels')

% %{
% erEnd = earlyRefsCutoff_samp+fade_samp_ER2LR+directSoundLength_samp-fade_samp_DS2ER;

subplot(2,4,7)
% figure;
plot(t_interp,squeeze(mean(rms(srirs_interp(erEnd:end,:,:)),2)));
hold on
plot(t_input,squeeze(mean(rms(srirs_input(erEnd:end,:,:)),2)));
title('Late Reverb RMS All channels')



subplot(2,4,8)
% figure;
plot(t_interp,mean(squeeze(rms(abs(squeeze(srirs_interp(1:directSoundLength_samp,:,:))))./rms(abs(squeeze(srirs_interp(directSoundLength_samp:end,:,:)))))));
hold on;plot(t_input,mean(squeeze(rms(abs(squeeze(srirs_input(1:directSoundLength_samp,:,:))))./rms(abs(squeeze(srirs_input(directSoundLength_samp:end,:,:)))))));
title('DRR All channels')
%}

%%

%% SPECTROGRAM
figure

irToPlot = 3;
channelToPlot = 24;
subplot(2,2,1)
spectrogram(srirs_interp(1:0.3*fs,channelToPlot,irToPlot),kaiser(256,5),220/2,512,fs,'yaxis');
title(['interpolated,ir=',num2str(irToPlot),', plotted channel = ',num2str(channelToPlot)]);
subplot(2,2,3)
% spectrogram(srirs_input(1:0.3*fs,channelToPlot,irToPlot),kaiser(256,5),220/2,512,fs,'yaxis');
title('correct')
% subplot(3,1,3)
% spectrogram(srirs_input(1:0.6*fs,channelToPlot,irToPlot)-srirs_interp(1:0.6*fs,channelToPlot,irToPlot),kaiser(256,5),220/2,512,fs,'yaxis');

% freq plot

% figure;
subplot(2,2,2)
freqplot(srirs_interp(directSoundLength_samp:earlyRefsCutoff_samp,channelToPlot,irToPlot),fs);
hold on
% freqplot(srirs_input(directSoundLength_samp:earlyRefsCutoff_samp,channelToPlot,irToPlot),fs);
title(['interpolated, ir=',num2str(irToPlot),', plotted channel = ',num2str(channelToPlot)]);
legend('interp','correct')

irToPlot2 = 9;
channelToPlot2 = 24;

subplot(2,2,4)
hold off
freqplot(srirs_interp(directSoundLength_samp:earlyRefsCutoff_samp,channelToPlot2,irToPlot2),fs);
hold on 
% freqplot(srirs_input(directSoundLength_samp:earlyRefsCutoff_samp,channelToPlot2,irToPlot2),fs);
title(['interpolated, ir=',num2str(irToPlot2),', plotted channel = ',num2str(channelToPlot2)]);
legend('interp','correct')


%% time domain plots

%{
channelToPlot = 25; % good to check not just W channel (1) to see if the higher order channels are good too
plot_length = 300; % plot length in ms

plot_length_samples = plot_length / 1000 * fs;
[X,Y] = meshgrid(0:5:500,1/48:1/48:plot_length);
% irTrunc = srirs_input;
% figure;
% surf(X,Y,real(10*log10(abs(squeeze(irTrunc(1:plot_length_samples,channelToPlot,:))))),'EdgeColor','none');
% 
% ylabel('Time (ms)');
% c = flip(bone);
% caxis([-35 0])
% colormap(c);
% h = colorbar;
% ylabel(h, 'Amplitude (dB)')
% view(0,90);
% set(gcf, 'Color', 'w');
% pbaspect([1.7 1 1]);
% set(gca, 'fontsize', 12);
% title(['time-domain - original, plotted channel = ',num2str(channelToPlot)]);


% interpolated

irTrunc = srirs_interp;
% plot_length = 40; % plot length in ms
% plot_length_samples = plot_length / 1000 * fs;
% [X,Y] = meshgrid(0:5:500,1/48:1/48:plot_length);

figure;
surf(X,Y,real(10*log10(abs(squeeze(irTrunc(1:plot_length_samples,channelToPlot,:))))),'EdgeColor','none');

ylabel('Time (ms)');
c = flip(bone);
caxis([-35 0])
colormap(c);
h = colorbar;
ylabel(h, 'Amplitude (dB)')
view(0,90);
set(gcf, 'Color', 'w');
pbaspect([1.7 1 1]);
set(gca, 'fontsize', 12);
title(['time-domain - interpolated, plotted channel = ',num2str(channelToPlot)]);




%}







%% DOA
%{
% % original
% irTrunc = srirs_input;
% 
% Tuneable parameters:
degreeResolution = 5;
order = 4;
nSrc = 7;
% numSamps = 48; % debugging
numSamps = 4800*5; % 0.5 seconds
highPassFilterFreq = 3000;
% 
        kappa = 40;
% 
% 
% % aziElev2aziPolar = @(dirs) [dirs(:,1) pi/2-dirs(:,2)]; % function to convert from azimuth-inclination to azimuth-elevation
grid_dirs = grid2dirs(degreeResolution,degreeResolution,0,0); % Grid of directions to evaluate DoA estimation
P_src = diag(ones(numSamps,1));
P_diff = 1; % unit power for the diffuse sound
[~,filtHi,~] = ambisonic_crossover(highPassFilterFreq,fs);
% 
% doa_est = zeros(nSrc,2,length(irTrunc(1,1,:)));
% doa_est_P = zeros(nSrc,length(irTrunc(1,1,:)));
% 
% for i = 1:length(irTrunc(1,1,:))
%     Y_src = filter(filtHi,1,irTrunc(1:numSamps,:,i)); % high pass filter
%     
%     stVec = Y_src';
%     sphCOV = stVec*P_src*stVec' + P_diff*eye((order+1)^2)/(4*pi);
%     
%     % DoA estimation
%     [~, est_dirs_pwd,est_dirs_P] = sphPWDmap(sphCOV, grid_dirs, nSrc,kappa);
%     
%     % convert to degs from rads
%     est_dirs_pwd = est_dirs_pwd*180/pi;
%     
%     % flip -ve values near -180 to +ve values
%     negativeFlipLimit = -170;
%     for j = 1:length(est_dirs_pwd(:,1))
%         for k = 1:length(est_dirs_pwd(1,:))
%             if est_dirs_pwd(j,k) < negativeFlipLimit
%                 est_dirs_pwd(j,k) = est_dirs_pwd(j,k) + 360;
%             end
%         end
%     end
%     doa_est(:,:,i) = est_dirs_pwd;
%     doa_est_P(:,i) = est_dirs_P;
% end
% 
% normalized_doa_est_P = doa_est_P ./ max(doa_est_P,[],1);
% normalized_doa_est_P_dB = mag2db(normalized_doa_est_P)/2;
% 
% plot_thresh = -10;
% normalized_doa_est_P_dB( normalized_doa_est_P_dB < plot_thresh ) = plot_thresh;
% 
% c_truncation = 20;
% cmap = flip(parula(256));
% c = cmap(c_truncation:end,:); % truncate yellow
% 
% doa_01 = rescale(normalized_doa_est_P_dB, 'InputMin',plot_thresh);
% 
% cspace = linspace(0,1,size(c,1));
% doa_color(:,:,1) = interp1(cspace, c(:,1), doa_01);
% doa_color(:,:,2) = interp1(cspace, c(:,2), doa_01);
% doa_color(:,:,3) = interp1(cspace, c(:,3), doa_01);
% 
%     figure;
% for i = nSrc:-1:1
%     s = scatter(0:5:500,squeeze(doa_est(i,1,:)),25,...
%         squeeze(doa_color(i,:,:)),'filled');
%     hold on
% end
% 
% ylabel('Azimuth (°)');
% 
% xlabel('<-- Room        Measurement position (cm)        Hallway -->');
% 
% ylim([negativeFlipLimit (negativeFlipLimit+360)]);yticks(-180:45:180);
% xlim([0 500]);xticks(0:50:500);
% colormap(c);
% 
% k = colorbar;
% ylabel(k,'Normalised power (dB)');
% caxis([plot_thresh 0]);
% set(gcf, 'Color', 'w');
% pbaspect([1.7 1 1]);
% box on;grid on;
% 
% title('DoA - original');
% 
% 
% 










%% Interpolated DOA

irTrunc = srirs_interp;



% degreeResolution = 5;
% order = 4;
% nSrc = 7;
% % numSamps = 48; % debugging
% numSamps = 4800*5; % 0.5 seconds
% highPassFilterFreq = 3000;


doa_est = zeros(nSrc,2,length(irTrunc(1,1,:)));
doa_est_P = zeros(nSrc,length(irTrunc(1,1,:)));

for i = 1:length(irTrunc(1,1,:))
    Y_src = filter(filtHi,1,irTrunc(1:numSamps,:,i)); % high pass filter
    
    stVec = Y_src';
    sphCOV = stVec*P_src*stVec' + P_diff*eye((order+1)^2)/(4*pi);
    
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

% cmap = flip(parula(256));
% c = cmap(c_truncation:end,:); % truncate yellow

doa_01 = rescale(normalized_doa_est_P_dB, 'InputMin',plot_thresh);

cspace = linspace(0,1,size(c,1));
doa_color(:,:,1) = interp1(cspace, c(:,1), doa_01);
doa_color(:,:,2) = interp1(cspace, c(:,2), doa_01);
doa_color(:,:,3) = interp1(cspace, c(:,3), doa_01);

    figure;
for i = nSrc:-1:1
    s = scatter(0:5:500,squeeze(doa_est(i,1,:)),25,...
        squeeze(doa_color(i,:,:)),'filled');
    hold on
end

ylabel('Azimuth (°)');

xlabel('<-- Room        Measurement position (cm)        Hallway -->');

ylim([negativeFlipLimit (negativeFlipLimit+360)]);yticks(-180:45:180);
xlim([0 500]);xticks(0:50:500);
colormap(c);

k = colorbar;
ylabel(k,'Normalised power (dB)');
caxis([plot_thresh 0]);
set(gcf, 'Color', 'w');
pbaspect([1.7 1 1]);
box on;grid on;

title(['DoA - interpolated. ']);%Measurement reduction = ',num2str(measReduction),'. interpType = ',interpType]);


%}

%% Determine the direct sound direction --- the points to be interpolated
figure
m = SRIR_Measurement();
m.setShMicArray(4, 0.032);
m.fs = fs;

numMeasurements = size(srirs_input,3);


directSoundDirections = zeros(numMeasurements, 3);
aziMeasuredDeg = zeros(numMeasurements, 1);
for iMeas = 1:numMeasurements
    m.srir = srirs_input(1:directSoundLength_samp, :, iMeas);
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

% Plot the DoAs
% figure
subplot(1,2,1)
quiver3(pos_input(:, 1) / 100, pos_input(:, 2) / 100, pos_input(:, 3) / 100, ...
    directSoundDirections(:, 1), directSoundDirections(:, 2) ,directSoundDirections(:, 3),0)

% text(pos_input(:, 1)/ 100, pos_input(:, 2)/ 100 ,pos_input(:, 3)/ 100, num2str((1:numMeasurements)'))
% hold on
% scatter3(loudspeakerPosition_cm(1) / 100, loudspeakerPosition_cm(2) / 100 , ...
%     loudspeakerPosition_cm(3) / 100, 20, 'rx')
axis equal
% zlim([1 2]);xlim([-.5 5.5]);ylim([-1.5 1.5]);
xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
title('Direct Sound DOAs - Positions to be Interpolated')
legend('Direct Sound DoA at the receiver', 'location', 'south')


%%%%%%%%%%%%%%%%%%% Compute DoA after Interpolation for checking
mInterp = SRIR_Measurement();
mInterp.fs = fs;
mInterp.setShMicArray(4, 0.032);
opts.analysisMethod = 'iv';
numInterpolatedPoints = size(srirs_interp,3);

directSoundDirectionsInterp = zeros(numInterpolatedPoints, 3);
for iInterp = 1:numInterpolatedPoints
    mInterp.srir = srirs_interp(1:directSoundLength_samp, :, iInterp);%hAllDirectSoundInterpolated
    a = SRIR_Analysis(mInterp, opts);
    a.run();
    
    idxNonZero = find(sum(a.doa~=0, 2));
    
    doaInterp = mean(a.doa(idxNonZero, :));
    doaInterp = doaInterp ./ (vecnorm(doaInterp, 2, 2));
    directSoundDirectionsInterp(iInterp, :) = doaInterp;
end

% Plot the DoA after interpolatioon
% figure
subplot(1,2,2)
quiver3(pos_interp(:, 1) / 100, pos_interp(:, 2) / 100, pos_interp(:, 3) / 100, ...
    directSoundDirectionsInterp(:, 1), directSoundDirectionsInterp(:, 2) ,directSoundDirectionsInterp(:, 3),0)
axis equal
% zlim([1 2]);xlim([-.5 5.5]);ylim([-1.5 1.5]);
% hold on
% scatter3(loudspeakerPosition_cm(1) / 100, loudspeakerPosition_cm(2) / 100 , ...
%     loudspeakerPosition_cm(3) / 100, 20, 'rx')
xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
title(['Direct Sound DOAs - Interpolated versions.'])
legend('Direct Sound DoA at the receiver', 'location', 'south')

