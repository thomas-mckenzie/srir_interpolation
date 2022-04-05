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

% choose new inter-measurement distance 
resolution_new = 5;

% minPhase, rotationOnly or fixedSpectrum or meanSpectrum (direct sound
% interpolation method)
INTERPOLATION_MODE_DS = 'meanSpectrum';

% directSoundLength_samp = 200; % empirically chosen
% fade_samp_DS2ER = 20; % DS direct sound ER early refs LR late reverb
% fade_samp_ER2LR = sofa.Data.SamplingRate / 100; % 480 samples
% 
% % SH order for interpolated signals
% NshInterp = 4;

fs = sofa.Data.SamplingRate;

%% REFERENCE DATASET (for evaluating the method)
% CORRECT = non-interpolated full resolution (for evaluating the method)
hCorrect = permute(sofa.Data.IR,[3,2,1]);
positionsCorrect = round(sofa.ListenerPosition*100); % in cm not m

%% MAKE SPARSE DATASET (for evaluating the method) 
% The amount to which to downsample the original measurements
% (2 would mean using 1 of every 2, 10 means 1 in every 10,
% 25 means 1 in every 25, so for 101 measurements it would just be using 5)
measReduction = 5;

% get the RIR dataset (for evaluating the method) - this would normally
% just be: 
% srirs_input = permute(sofa.Data.IR,[3,2,1]);
% srirs_input = hCorrect(:,:,1:measReduction:end);
% 
% % reduce the reference dataset (for evaluating the method)  --- this should
% % probably normally just be:
% % pos_input = round(sofa.ListenerPosition*100); 
% pos_input = positionsCorrect(1:measReduction:end,:);
% 

srirs_input = hCorrect(:,:,1:measReduction:end);

% reduce the reference dataset (for evaluating the method)  --- this should
% probably normally just be:
% pos_input = round(sofa.ListenerPosition*100); 
pos_input = positionsCorrect(1:measReduction:end,:);


%% RUN THE FUNCTION

[srirs_interp,pos_interp] = interpolate_SRIRs(srirs_input,pos_input,resolution_new,fs,INTERPOLATION_MODE_DS);


%% evaluate
interpType = 'sectors, no time windows';

run analyse_interpolated_SRIRs_function_noCorrect.m

%% Save to sofa
%{
SOFAstart()

Obj = sofa;
Obj.Data.IR = permute(srirs_interp, [3 2 1]);

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
switch INTERPOLATION_MODE_DS
    case 'meanSpectrum'
        SOFAfn = fullfile(['srirInterp_ms_',num2str(measReduction),'.sofa']);
    case 'minPhase'
        SOFAfn = fullfile(['srirInterp_mp_',num2str(measReduction),'.sofa']);
end

disp(['Saving:  ' SOFAfn]);
Obj = SOFAsave(SOFAfn, Obj, 1);
%}
% end