% Analysis script for interpolate_SRIRs_source.m
%
% The user is directed to the following paper:
% McKenzie, T. and Schlecht, S. J. "Source position interpolation of 
%     spatial room impulse responses" AES 154th Convention, Espoo, Helsinki, 
%     Finland, 2023. 
%
% Thomas McKenzie, University of Edinburgh, 2023. thomas.mckenzie@ed.ac.uk


%% Time-domain plot
t_interp = 1:size(srirs_interp,3);
figure;
plot(t_interp,db(rms((squeeze(srirs_interp(:,1,:,1))))),'LineWidth',2,'Marker','^','MarkerSize',10,'LineStyle','none');
hold on
plot(t_interp,db(rms((squeeze(srirs_interp(:,1,:,end))))),'LineWidth',2,'Marker','^','MarkerSize',10,'LineStyle','none');
numInterpResp = size(srirs_interp,4)-2;
for i = 2:size(srirs_interp,4)-1
    plot(t_interp,db(rms((squeeze(srirs_interp(:,1,:,i))))),':','Color',[((i-1)/numInterpResp - 1/numInterpResp/2) 0 1-((i-1)/numInterpResp-1/numInterpResp/2)],'Marker','^','MarkerSize',10,'LineStyle','none');
end

xlim([0.5 size(srirs_interp,3)+0.5])
ylabel('Magnitude (dB)'); xlabel('Receiver Number')
pbaspect([2.2 1 1])
ylim([-49 -38]);box on
title(['RMS W channel: ', methodType, ' interpolation method'])

%% Frequency-domain response
irToPlot = 1; channelToPlot = 1;
figure;
freqplot_smooth(srirs_interp(:,channelToPlot,irToPlot,1),fs,2,'-', 2);
hold on
freqplot_smooth(srirs_interp(:,channelToPlot,irToPlot,end),fs,2,'-', 2);
for i = 2:size(srirs_interp,4)-1
    freqplot_smooth(srirs_interp(:,channelToPlot,irToPlot,i),fs,2,'k:',1);
end

pbaspect([2.2 1 1]); ylim([-15 25]);
title(['Frequency response W channel: ', methodType, ' interpolation method'])

%% Direct sound DoA

figure;
hold on;grid on;

m = SRIR_Measurement(); m.setShMicArray(4, 0.032); m.fs = fs;
numMeasurements = size(srirs_interp,3);
directSoundLength_samp = 150;

% DoA, first source position
directSoundDirections = zeros(numMeasurements, 3);
aziMeasuredDeg = zeros(numMeasurements, 1);
for iMeas = 1:numMeasurements
    m.srir = squeeze(srirs_interp(1:directSoundLength_samp, :, iMeas,1));
    opts.analysisMethod = 'iv';
    a = SRIR_Analysis(m, opts);    a.run();

    idxNonZero = find(sum(a.doa~=0, 2));
    doa = mean(a.doa(idxNonZero, :));
    doa = doa ./ vecnorm(doa, 2, 2);

    directSoundDirections(iMeas, :) = doa;
    aziMeasuredDeg(iMeas) = atan2(doa(2), doa(1)) * 180 / pi;
    if aziMeasuredDeg(iMeas)<0
        aziMeasuredDeg(iMeas) = aziMeasuredDeg(iMeas) + 360;
    end
end
quiver3(recPosns(:, 1) / 100, recPosns(:, 2) / 100, recPosns(:, 3) / 100, ...
    directSoundDirections(:, 1), directSoundDirections(:, 2) ,directSoundDirections(:, 3),0,'Linewidth',2)

% DoA, second source position
directSoundDirections = zeros(numMeasurements, 3);
aziMeasuredDeg = zeros(numMeasurements, 1);
for iMeas = 1:numMeasurements
    m.srir = squeeze(srirs_interp(1:directSoundLength_samp, :, iMeas,end));
    a = SRIR_Analysis(m, opts);    a.run();

    idxNonZero = find(sum(a.doa~=0, 2));
    doa = mean(a.doa(idxNonZero, :));
    doa = doa ./ vecnorm(doa, 2, 2);

    directSoundDirections(iMeas, :) = doa;
    aziMeasuredDeg(iMeas) = atan2(doa(2), doa(1)) * 180 / pi;
    if aziMeasuredDeg(iMeas)<0
        aziMeasuredDeg(iMeas) = aziMeasuredDeg(iMeas) + 360;
    end
end
quiver3(recPosns(:, 1) / 100, recPosns(:, 2) / 100, recPosns(:, 3) / 100, ...
    directSoundDirections(:, 1), directSoundDirections(:, 2) ,directSoundDirections(:, 3),0,'Linewidth',2)

% DoA, interpolated source positions
for iIntSrc = 2:size(srirs_interp,4)-1
    directSoundDirections = zeros(numMeasurements, 3);
    aziMeasuredDeg = zeros(numMeasurements, 1);
    for iMeas = 1:numMeasurements
        m.srir = squeeze(srirs_interp(1:directSoundLength_samp, :, iMeas,iIntSrc));
        a = SRIR_Analysis(m, opts);        a.run();

        idxNonZero = find(sum(a.doa~=0, 2));
        doa = mean(a.doa(idxNonZero, :));
        doa = doa ./ vecnorm(doa, 2, 2);

        directSoundDirections(iMeas, :) = doa;
        aziMeasuredDeg(iMeas) = atan2(doa(2), doa(1)) * 180 / pi;
        if aziMeasuredDeg(iMeas)<0
            aziMeasuredDeg(iMeas) = aziMeasuredDeg(iMeas) + 360;
        end
    end
    quiver3(recPosns(:, 1) / 100, recPosns(:, 2) / 100, recPosns(:, 3) / 100, ...
        directSoundDirections(:, 1), directSoundDirections(:, 2) ,directSoundDirections(:, 3),0,'Color',[((iIntSrc-1)/numInterpResp - 1/numInterpResp/2) 0 1-((iIntSrc-1)/numInterpResp-1/numInterpResp/2)])
    axis equal
    xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
end

% plot actual source positions (and interpolated 'correct' source positions)
scatter3(sofa1.SourcePosition(src1,1),sofa1.SourcePosition(src1,2),sofa1.SourcePosition(src1,3),50,[0, 0.4470, 0.7410],'LineWidth',2,'Marker','^')
scatter3(sofa1.SourcePosition(src2,1),sofa1.SourcePosition(src2,2),sofa1.SourcePosition(src2,3),50,[0.8500, 0.3250, 0.0980],'LineWidth',2,'Marker','^')
for i = 2:size(srirs_interp,4)-1
    scatter3(pos_interp(i,1)/100,pos_interp(i,2)/100,pos_interp(i,3)/100,50,[((i-1)/numInterpResp - 1/numInterpResp/2) 0 1-((i-1)/numInterpResp-1/numInterpResp/2)],'Marker','^');
end

xlim([1 7.87]);ylim([1 5]);zlim([1 2]);view([0 90]);box on;
title(['Direct sound DoA: ', methodType, ' interpolation method'])
