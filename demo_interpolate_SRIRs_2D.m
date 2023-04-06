% The user is directed to the ICA 2022 paper for details:
% ﻿McKenzie, T., Meyer-Kahlen, N., Daugintis, R., McCormack, L., Schlecht, S. 
% J., & Pulkki, V. (2022). Perceptual interpolation and rendering of coupled 
% room spatial room impulse responses. International Congress on Acoustics, 
% Korea. 
% 
% This demo script uses a 2D dataset of SRIRs.  
% 
% Thomas McKenzie, 2022. thomas.mckenzie@aalto.fi / tom.mckenzie07@gmail.com

% positions are in the motive system with the "z up" setting which is the 
% same, right handed, x in the front system that sofa and the 6DoFconv VST
% use x front y left z up

% close all
clear
clc

sofa = SOFAload('arni_dataset_no_interp.sofa');

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

[srirs_interp,pos_interp] = interpolate_SRIRs(srirs_input,pos_input,resolution_new,fs,INTERPOLATION_MODE_DS);

%% evaluate
interpType = 'sectors, no time windows';

run analyse_interpolated_SRIRs_function_noCorrect.m


%% Save to sofa
% %{
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
% switch INTERPOLATION_MODE_DS
%     case 'meanSpectrum'
        SOFAfn = fullfile(['srirInterp_2d_noWin.sofa']);
%     case 'minPhase'
%         SOFAfn = fullfile(['srirInterp_mp_',num2str(measReduction),'.sofa']);
% end

disp(['Saving:  ' SOFAfn]);
Obj = SOFAsave(SOFAfn, Obj, 1);
%}
% end
