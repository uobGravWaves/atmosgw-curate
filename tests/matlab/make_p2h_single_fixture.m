%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB reference output for p2h with single-precision input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Writes ../data/matlab_p2h_single.mat: single-precision pressures (hPa) and
%p2h's heights for them (km, double), for tests/test_atmosphere.py. With single
%input, p2h computes in single precision (as for MLS, whose pressures are
%single).
%
%How to run (set the two variables first):
%   MatlabFunctionsDir    = '<path to a checkout of corwin365/MatlabFunctions>';
%   MatlabFunctionsCommit = '<its commit hash>';
%   run('<full path>/make_p2h_single_fixture.m')
%It refuses to overwrite an existing output file.
%
%change log:
%  2026/11/24 [Claude] created
%  2026/11/24 [Claude] machine name no longer stored in the output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('MatlabFunctionsDir','var') || ~exist('MatlabFunctionsCommit','var')
  error('Set MatlabFunctionsDir and MatlabFunctionsCommit before running (see header)');
end
ThisDir = fileparts(mfilename('fullpath'));
OutFile = fullfile(ThisDir,'..','data','matlab_p2h_single.mat');
if exist(OutFile,'file'); error('%s already exists: delete it first if you mean to remake it',OutFile); end
addpath(fullfile(MatlabFunctionsDir,'GeophysicalAndTime'))
rng(20261125,'twister')

%random pressures from 1e-3 to 1000 hPa, and a grid of 12 levels per decade from
%1000 hPa (like the MLS temperature levels)
P = single([10.^(6.*rand(1,300) - 3), 1000.*10.^(-(0:54)./12)]);
Z = p2h(P);

Meta = struct('Created',datestr(now,'yyyy-mm-dd HH:MM:SS'),'Matlab',version, ...
              'MatlabFunctionsCommit',MatlabFunctionsCommit,'Seed',20261125);
save(OutFile,'Meta','P','Z','-v7');
fprintf('saved %s (class P %s, class Z %s)\n',OutFile,class(P),class(Z));
