%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make MATLAB reference outputs for atmosgw.curate.analysis.detrend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Writes OutFile (<reference_dir>/detrend_sgolay_ref.mat): what the group's
%MATLAB detrending does to real limb profiles and a small random array,
%compared against by tests/test_detrend_matlab.py. Inputs are saved with the
%outputs, so the Python test needs no other file.
%
%Covered:
%  group's own code:  func_filter_sgolay, func_filter_hindley23 (private
%                     functions of gwanalyse_limb, corwin365/MatlabFunctions),
%                     called exactly as gwanalyse_limb calls them (Data.Tp = Data.Temp)
%  third party:       smoothn (A. Gannepali, File Exchange 725), box and gaussian
%  (the cross-track polynomial is compared against prep_airs_3d_reference.mat)
%
%Inputs: the profiles saved in gwanalyse_limb_reference.mat (HIRDLS and SABER
%for 2007-01-01 from get_limbsounders; SABER_origz on its original heights;
%HIRDLS_PW the Hindley23 planetary-wave files).
%
%How to run (set the variables first; the run records them in the file):
%   MatlabFunctionsDir = '<path to checkouts of corwin365/MatlabFunctions-Common and -Personal>';
%   CommonCommit       = '<commit of Common>';   PersonalCommit = '<commit of Personal>';
%   LimbReference      = '<full path>/gwanalyse_limb_reference.mat';
%   OutFile            = '<full path>/detrend_sgolay_ref.mat';
%   run('<full path>/make_detrend_reference.m')
%It refuses to overwrite an existing output file.
%
%change log:
%  2026/11/28 [Claude] created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Needed = {'MatlabFunctionsDir','CommonCommit','PersonalCommit','LimbReference','OutFile'};
for iN=1:1:numel(Needed)
  if ~exist(Needed{iN},'var'); error('Set %s before running (see header)',Needed{iN}); end
end
if exist(OutFile,'file'); error('%s already exists: delete it first if you mean to remake it',OutFile); end
addpath(genpath(MatlabFunctionsDir))
%the filters are private functions of gwanalyse_limb: put their folder on the path
addpath(fullfile(fileparts(which('gwanalyse_limb')),'gwanalyse_limb'))

Limb = load(LimbReference,'In');
Ref  = struct();

%% Savitzky-Golay: [Instrument, SGLength (km), SGOrder]
Cases = {'HIRDLS',18,2; 'SABER',18,2; 'SABER_origz',18,2; 'HIRDLS',5,2; 'HIRDLS',10,3};
for iC=1:1:size(Cases,1)
  Data = struct('Temp',Limb.In.(Cases{iC,1}).Temp,'Alt',Limb.In.(Cases{iC,1}).Alt);
  Data.Tp = Data.Temp;                                         %as gwanalyse_limb
  Settings = struct('SGLength',Cases{iC,2},'SGOrder',Cases{iC,3},'Verbose',0);
  Out  = func_filter_sgolay(Data,Settings);
  Name = sprintf('SG_%s_%d_%d',Cases{iC,1},Cases{iC,2},Cases{iC,3});
  Ref.(Name) = struct('Temp',Data.Temp,'Alt',Data.Alt,'Tp',Out.Tp,'BG',Out.BG, ...
                      'SGLength',Cases{iC,2},'SGOrder',Cases{iC,3});
end

%% Hindley23: copy of the stored residual and planetary waves, with and without outlier removal
In = Limb.In.HIRDLS_PW;
for OutRem=[false,true]
  Data = In; Data.Tp = [];
  [Out,PW] = func_filter_hindley23(Data,struct('H23_OutRem',OutRem));
  Name = sprintf('H23_OutRem%d',OutRem);
  Ref.(Name) = struct('Tp',Out.Tp,'BG',Out.BG);
end
Ref.H23_In = struct('Temp_Residual',In.Temp_Residual,'Temp_PW',In.Temp_PW,'Alt',In.Alt,'Lat',In.Lat,'Lon',In.Lon);

%% smoothn: box and gaussian, [cross x along x level]-like random array, with and without a NaN
rng(1); X = 250 + 10.*randn(14,11,7);
XNaN = X; XNaN(3,4,2) = NaN;
Ref.Smooth.X          = X;
Ref.Smooth.XNaN       = XNaN;
Ref.Smooth.Box_3_5_1  = smoothn(X,[3,5,1]);
Ref.Smooth.Box_3_3_3  = smoothn(X,[3,3,3]);
Ref.Smooth.BoxNaN     = smoothn(XNaN,[3,5,1]);
Ref.Smooth.Gauss_5_3_3     = smoothn(X,[5,3,3],'gaussian');               %default std 0.65
Ref.Smooth.Gauss_5_5_1_std = smoothn(X,[5,5,1],'gaussian',[1,0.5,0.8]);
Ref.Smooth.Gauss_std       = [1,0.5,0.8];

%% record how this was made
Ref.Meta.Created        = datestr(now,'yyyy-mm-dd HH:MM:SS');
Ref.Meta.Matlab         = version;
Ref.Meta.CommonCommit   = CommonCommit;
Ref.Meta.PersonalCommit = PersonalCommit;
Ref.Meta.LimbReference  = LimbReference;

save(OutFile,'Ref','-v7');
fprintf('Wrote %s\n',OutFile);
