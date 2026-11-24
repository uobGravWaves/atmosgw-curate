%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write the synthetic input files used by the MATLAB get_limbsounders tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%The MATLAB test (MatlabFunctions/Tests/test_get_limbsounders.m) writes small
%synthetic input files to a temporary folder, runs get_limbsounders on them,
%saves the outputs in get_limbsounders_reference.mat, and deletes the files.
%To compare the Python port with those saved outputs, it must read the very
%same input files. This script writes them, with the same code as the test's
%make_synthetic_files (copied unchanged below), to a folder that is kept.
%
%How to run (set the variable first):
%   SynDir = '<reference folder>/get_limbsounders_synthetic';
%   run('<full path>/make_limb_synthetic_inputs.m')
%It refuses to write into a folder that already exists.
%
%The files: 2007-01-01 and 2007-01-02, 8 profiles x 21 levels per day.
%Profile k has heights 10:2:50 km + 0.3k, Temp = 200 + 1.5*Alt, Lat
%linspace(-60,60,8), Lon linspace(-157.5,157.5,8), time (k-0.5)/8 of the day.
%Variants: synlin (plain), synbad (outliers), syn3d (3-D extra variable),
%synpw (PW format), synext (extra variables of different shapes), synlon
%(longitudes outside [-180,180]) and sampler/ (sampler format, day 1 only).
%
%change log:
%  2026/11/24 [Claude] created (make_synthetic_files copied from
%                      test_get_limbsounders.m, 2026/11/22 version)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('SynDir','var'); error('Set SynDir before running (see header)'); end
if exist(SynDir,'dir'); error('%s already exists: choose a new folder',SynDir); end
mkdir(SynDir);
make_synthetic_files(SynDir);
fprintf('wrote synthetic files to %s\n',SynDir);


function make_synthetic_files(SynDir)
%write small synthetic input files with known answers, for 2007-01-01 and 2007-01-02.
%8 profiles x 21 levels per day. Profile k: heights 10:2:50 km + 0.3k (so interpolation is
%non-trivial), Temp = 200 + 1.5*Alt (linear, so interpolation is exact), Lat linspace(-60,60,8),
%Lon linspace(-157.5,157.5,8), time (k-0.5)/8 of the day
d0 = datenum(2007,1,1); NP = 8; Lev = 10:2:50; NL = numel(Lev);
mkdir(fullfile(SynDir,'sampler'));
for iDay=0:1
  Tag = sprintf('2007d%03d',iDay+1);
  S = struct();
  S.Alt  = 0.3.*(1:NP)' + Lev;                       %[NP x NL] (implicit expansion)
  S.Lat  = repmat(linspace(-60,60,NP)',1,NL);
  S.Lon  = repmat(linspace(-157.5,157.5,NP)',1,NL);
  S.Time = repmat(d0+iDay+((1:NP)'-0.5)./NP,1,NL);
  S.Temp = 200 + 1.5.*S.Alt;
  save(fullfile(SynDir,['synlin_',Tag,'.mat']),'-struct','S');

  %outliers on day 1
  B = S; if iDay == 0; B.Lat(3,5) = 95; B.Temp(6,7) = 50; end
  save(fullfile(SynDir,['synbad_',Tag,'.mat']),'-struct','B');

  %3-D additional variable
  T = S; T.Extra3 = S.Temp.*reshape(1:3,1,1,3);
  save(fullfile(SynDir,['syn3d_',Tag,'.mat']),'-struct','T');

  %PW-format file: Temp_PW(:,:,k) = k*Alt; outlier latitude on day 1
  P = rmfield(S,'Temp'); P.Temp_Residual = 0.1.*S.Alt; P.Temp_PW = S.Alt.*reshape(1:4,1,1,4);
  if iDay == 0; P.Lat(3,5) = 95; end
  P.Meta = 7;  %scalar metadata, only loaded if requested
  save(fullfile(SynDir,['synpw_',Tag,'.mat']),'-struct','P');

  %extra variables of different shapes
  %(PerLev is non-integer so that it is interpolated as a physical variable; FlagVar is an integer flag 0,1,2 by level)
  X = S; X.PerPointT = S.Temp'; X.PerProf = 10.*(1:NP)'; X.PerLev = Lev + 0.5; X.Meta = 7; X.TempMeta = 3;
  X.FlagVar = repmat(mod(0:NL-1,3),NP,1);
  save(fullfile(SynDir,['synext_',Tag,'.mat']),'-struct','X');

  %longitudes outside [-180,180]: profiles 1-4 at 190, -190, 350, -350 (expected -170, 170, -10, 10)
  W = S; W.Lon = repmat([190;-190;350;-350;0;10;20;30],1,NL);
  save(fullfile(SynDir,['synlon_',Tag,'.mat']),'-struct','W');

  %sampler format: struct Sampled_Data, fields stored [levels x profiles]
  if iDay == 0
    Sampled_Data = struct('Alt',S.Alt','Lat',S.Lat','Lon',S.Lon','Time',S.Time','T',S.Temp','Prs',h2p(S.Alt)','SourceProf',repmat(1:NP,NL,1));
    save(fullfile(SynDir,'sampler',['sampled_',num2str(d0),'.mat']),'Sampled_Data');
  end
end
end
