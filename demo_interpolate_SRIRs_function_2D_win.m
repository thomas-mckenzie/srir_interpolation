% close all
clear
clc

%%%%%%%%%%%%%%% 2D ----- NILS ARNI DATASET

% positions are in the motive system with the "z up" setting
% which is the same, right handed, x in the front system that
% sofa and the tv convolver use
% x front y left z up

path = 'C:\Users\nils\OneDrive - Aalto University\2021_transfer_plausibility\sofaFiles\nointerp\'
filename = 'm5_ls1.sofa'
sofa = SOFAload([path filename]);

%% configure

% choose new inter-measurement distance 
resolution_new = 10;

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

% get the RIR dataset
srirs_input = permute(sofa.Data.IR,[3,2,1]);
% get positions
pos_input = round(sofa.ListenerPosition*100); 

%% RUN THE FUNCTION

% response is partioned into regions containing proportions of energy.

region_limits_db = -5:-10:-50;

[srirs_interp,pos_interp] = ...
    interpolate_SRIRs_win(srirs_input,pos_input,resolution_new, ...
    fs,INTERPOLATION_MODE_DS, region_limits_db);

%% evaluate
%interpType = 'sectors, time windows';
%run analyse_interpolated_SRIRs_function_noCorrect.m

%% Save to sofa

SOFAstart()

Obj = sofa;
Obj.Data.IR = permute(srirs_interp, [3 2 1]);

Obj.ListenerPosition = pos_interp;
Obj.SourcePosition = [0, 0, 0];

% Update dimensions
Obj=SOFAupdateDimensions(Obj);

% Fill with attributes
Obj.GLOBAL_ListenerShortName = 'EM';
Obj.GLOBAL_History = ['created on ' string(datetime)];
Obj.GLOBAL_DatabaseName = 'none';
Obj.GLOBAL_ApplicationName = 'SOFA API';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
Obj.GLOBAL_Organization = 'Aalto Acoustics Lab';
Obj.GLOBAL_AuthorContact = 'thomas.mckenzie@aalto.fi';
Obj.GLOBAL_Comment = ' ';
Obj.GLOABL_Title =  'Responses on a grid, interpolated';

 SOFAfn = fullfile(['srirInterp_' filename]);

disp(['Saving:  ' SOFAfn]);
Obj = SOFAsave(SOFAfn, Obj, 1);
