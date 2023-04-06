% Demo script for interpolate_SRIRs_source.m 
% 
% The user is directed to the following paper: 
% McKenzie, T. and Schlecht, S. J. "Source position interpolation of 
%     spatial room impulse responses" AES 154th Convention, Espoo, Helsinki, 
%     Finland, 2023. 
% 
% This script uses the following dataset of spatial room impulse response
% measurements: 
% McKenzie, T., McCormack, L., and Hold, C., "Dataset of spatial room 
%     impulse responses in a variable acoustics room for six degrees-of-
%     freedom rendering and analysis," in arXiv preprint, pp. 1–3, 2021, 
%     doi:10.48550/arXiv.2111.11882.
% Dataset download link: https://doi.org/10.5281/zenodo.5720723
% 
% Thomas McKenzie, University of Edinburgh, 2023. thomas.mckenzie@ed.ac.uk

close all
clear
clc

% Load Variable Acoustics 6DoF Dataset (most reverberant set)
irPath1 = '/Users/tmckenzi/Documents/6DoF Dataset/6dof_SRIRs_eigenmike_SH/';
irName1 = '6DoF_SRIRs_eigenmike_SH_0percent_absorbers_enabled.sofa';

sofa1 = SOFAload([irPath1,irName1]);
fs = sofa1.Data.SamplingRate;
SOFAplotGeometry(sofa1);

%% Configure interpolation

division = 5; % divide into x amount of sections (ie 2: one new srir in the middle, 3: 2 new srirs, 4: 3 new srirs, 10: 9 new srirs)

% how many receiver positions (and which ones in index)
recPosnSel = 1:7; % all positions
% recPosnSel = [2,4,6,7]; % rectangle
% recPosnSel = 1:5; % line

% which two sources to use: in the 6DoF dataset, 1 is bottom right, 2 is bottom left, and 3 is upper middle
src1 = 1; src2 = 2;

%% Arrange Input DATASET (for evaluating the method)

src1_ind = recPosnSel*3-3+src1;
src2_ind = recPosnSel*3-3+src2;

sofaIRdata = permute(sofa1.Data.IR,[3,2,1]);

for i = 1:size(sofaIRdata,3) % pre-process of 6DoF SRIRs - 90 degree rotation
sofaIRdata(:,:,i) = rotateHOA_N3D(sofaIRdata(:,:,i),90,0,0);
end

h1 = sofaIRdata(:,:,src1_ind);
h2 = sofaIRdata(:,:,src2_ind);

positionsCorrect = round(sofa1.ListenerPosition(src1_ind,:)*100);

source1_position = round(sofa1.SourcePosition(src1_ind,:)*100);
source2_position = round(sofa1.SourcePosition(src2_ind,:)*100);

srirs_input1 = h1;
srirs_input2 = h2;

recPosns = positionsCorrect;
srirs_input = cat(4,srirs_input1, srirs_input2); % if interp all receiver posns, cat 
pos_input = [source1_position(1,:);source2_position(1,:)]; % assuming source position doesn't change with receiver position

%% RUN THE FUNCTION -- BASIC
[srirs_interp,pos_interp] = interpolate_SRIRs_source_basic(srirs_input,pos_input,division,fs);

%% Evaluate 
methodType = 'Basic';
run analyse_interpolate_SRIRs_source.m

%% RUN THE FUNCTION -- PROPOSED
[srirs_interp,pos_interp,dir_interp_src] = interpolate_SRIRs_source(srirs_input,pos_input,division,fs);

%% Evaluate 
methodType = 'Proposed';
run analyse_interpolate_SRIRs_source.m

%% Save to sofa

%{
% Need to decide what to save: ie to save either the rotating source for
% one receiver position, or a fixed (interpolated) source for moving receiver
% positions.
% fixed receiver position, moving source position:
receiver_position = 1;
srirs_interp = squeeze(srirs_interp(:,:,receiver_position,:));
% fixed source position (interpolated), moving receiver position:
% interp_source_position = round(size(srirs_interp,4)/2); % take middle one
% srirs_interp = squeeze(srirs_interp(:,:,:, interp_source_position)); 

SOFAstart()
Obj = sofa1;
Obj.Data.IR = permute(srirs_interp, [3 2 1]);

% Get listener and source position for the SOFA file:
% Correct VERSION - This is the actual right way to do it
% Obj.ListenerPosition = repmat(sofa1.ListenerPosition(inputsrir,:),division+1,1);
% Obj.SourcePosition = pos_interp/100;
% Compatibility VERSION - For compatibility with the Sparta 6DoFconv
plugin, which is only made for switching between listener positions
Obj.ListenerPosition = pos_interp/100;
Obj.SourcePosition = repmat(sofa1.ListenerPosition(inputsrir,:),division+1,1);

% useful for analysis - shows the movement of the direct path from source
% to receiver, but not for the SOFA file - this isn't the 'absolute' source
% view. --- paper shows the arrow plot from the interpolated source
% positions. 
% [dir_in_src_cart(1,:),dir_in_src_cart(2,:),dir_in_src_cart(3,:)] = sph2cart(dir_interp_src(1,:),dir_interp_src(2,:),1);

% Rotate source view
point1 = sofa1.SourceView; 
point2 = sofa2.SourceView;
[p1(1), p1(2) ] = cart2sph(point1(1),point1(2),point1(3));
[p2(1), p2(2) ] = cart2sph(point2(1),point2(2),point2(3));
p_diff = angdiff(p1,p2);
for i = 1:2
    if p_diff(i) == 0
      p_diff_div(:,i) = p1(i);
    else
    p_diff_div(:,i) = p1(i):-p_diff(i)/division:p2(i);
    end
end
[p3(:,1),p3(:,2),p3(:,3)] = sph2cart(p_diff_div(:,1),p_diff_div(:,2),1);

Obj.SourceView = p3;
Obj=SOFAupdateDimensions(Obj);

SOFAplotGeometry(Obj);

% Update sofa dimensions

% Fill with attributes
Obj.GLOBAL_ListenerShortName = 'EM';
Obj.GLOBAL_History = 'created on 07.04.2023';
Obj.GLOBAL_DatabaseName = 'none';
Obj.GLOBAL_ApplicationName = 'SOFA API';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
Obj.GLOBAL_Organization = 'Acoustics and Audio Group, University of Edinburgh';
Obj.GLOBAL_AuthorContact = 'thomas.mckenzie@ed.ac.uk';
Obj.GLOBAL_Comment = ' ';
Obj.GLOABL_Title =  'Source position interpolation';

% save the SOFA file
SOFAfn = fullfile(['srir_src_interp.sofa']);

disp(['Saving:  ' SOFAfn]);
Obj = SOFAsave(SOFAfn, Obj, 1);
%}


