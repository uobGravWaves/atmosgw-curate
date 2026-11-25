%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make MATLAB reference outputs for atmosgw.curate.swath and .solar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Writes ../data/matlab_swath_fixtures.mat: small synthetic swaths and what the
%MATLAB prep_airs_3d does to them, plus solar zenith angles and the AIRS
%day/night mask, compared against by tests/test_swath.py and
%tests/test_solar.py. Inputs are saved with the outputs.
%
%Covered:
%  Mapping Toolbox:   distance (spherical default: arc in degrees and azimuth)
%  group's own code:  prep_airs_3d regridding (InputStruct, no detrending), with
%                     and without vertical interpolation, by number and by
%                     spacing, linear and nearest; which_airs_retrieval
%                     (corwin365/MatlabFunctions)
%
%Swaths: three synthetic granules, built with reckon so their geometry is
%exact: mid-latitudes, one crossing the dateline, one near the AIRS turning
%latitude. 31 cross-track x 16 along-track points, 13 levels (AIRS-like:
%3 km spacing, then wider above 30 km).
%
%How to run (set the two variables first; the run records them in the file):
%   MatlabFunctionsDir    = '<path to a checkout of corwin365/MatlabFunctions>';
%   MatlabFunctionsCommit = '<its commit hash>';
%   run('<full path>/make_swath_fixtures.m')
%It refuses to overwrite an existing output file.
%
%change log:
%  2026/11/25 [Claude] created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('MatlabFunctionsDir','var') || ~exist('MatlabFunctionsCommit','var')
  error('Set MatlabFunctionsDir and MatlabFunctionsCommit before running (see header)');
end
ThisDir = fileparts(mfilename('fullpath'));
OutFile = fullfile(ThisDir,'..','data','matlab_swath_fixtures.mat');
if exist(OutFile,'file'); error('%s already exists: delete it first if you mean to remake it',OutFile); end
if ~license('test','map_toolbox'); error('Mapping Toolbox needed'); end
addpath(genpath(MatlabFunctionsDir))
rng(20261125,'twister')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% distance (arc length in degrees, and azimuth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DistA = [asind(2.*rand(200,1)-1), 360.*rand(200,1)-180];
DistB = [asind(2.*rand(200,1)-1), 360.*rand(200,1)-180];
%close pairs (swath-like spacing), the dateline, and identical points
Near = DistA(1:50,:) + 0.2.*randn(50,2); Near(:,1) = max(min(Near(:,1),90),-90);
DistA = [DistA; DistA(1:50,:); 10 179.9; -20 -179.95; 30 40];
DistB = [DistB; Near;          10 -179.9; -20 179.95; 30 40];
[DistArc,DistAz] = distance(DistA(:,1),DistA(:,2),DistB(:,1),DistB(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% synthetic swaths, regridded by prep_airs_3d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z  = [0 3 6 9 12 15 18 21 24 27 30 35 40];
NX = 31; NA = 16;
XT = sign((1:NX)-16).*14.*abs((1:NX)-16).^1.1;           %km from the middle, widening outwards
Geo = struct('Name',{'mid','dateline','turning'}, ...
             'Lat0',{40,-30,78},'Lon0',{10,178.8,-60},'Az0',{350,30,265});
Calls = struct('Name',{'linear','vertical','spacing','nearest'}, ...
               'Args',{{'VerticallyInterpolate',false}, ...
                       {}, ...
                       {'VerticallyInterpolate',false,'dXT',40,'dAT',30}, ...
                       {'VerticallyInterpolate',false,'Interpolant','nearest'}});
Swath = struct();
for iG = 1:numel(Geo)
  G = Geo(iG);
  %mid-track points every 18 km, heading changing slowly (a curved track)
  MidLat = NaN(1,NA); MidLon = NaN(1,NA); Head = G.Az0 + linspace(0,6,NA);
  MidLat(1) = G.Lat0; MidLon(1) = G.Lon0;
  for k = 2:NA; [MidLat(k),MidLon(k)] = reckon(MidLat(k-1),MidLon(k-1),km2deg(18),Head(k-1)); end
  [Lat,Lon] = reckon(repmat(MidLat,NX,1),repmat(MidLon,NX,1),km2deg(repmat(XT',1,NA)),repmat(Head+90,NX,1));
  T0 = (datenum(2010,6,1,3,0,0) - datenum(2000,1,1)).*86400;
  Time = T0 + repmat((0:NA-1).*2.6667,NX,1) + repmat((0:NX-1)'.*0.0222,1,NA);   %seconds since 2000-01-01
  [z,x,a] = ndgrid(Z,1:NX,1:NA);
  T = 200 + 0.8.*z + 3.*sin(2*pi*a/9).*cos(2*pi*x/17) + 0.02.*x.*a;
  T(6,20,10) = NaN;                                       %one missing value
  S = struct('l1_lat',Lat,'l1_lon',Lon,'l1_time',Time,'ret_temp',T,'ret_z',Z);
  Swath.(G.Name).In = S;
  for iC = 1:numel(Calls)
    [A,Sp,E,EI] = prep_airs_3d(0,0,'InputStruct',S,'NoDetrend',true,'NXT',16,'NAT',12,Calls(iC).Args{:});
    if E ~= 0; error('prep_airs_3d failed (%s, %s): %s',G.Name,Calls(iC).Name,EI); end
    Swath.(G.Name).(Calls(iC).Name) = struct('Lat',A.l1_lat,'Lon',A.l1_lon,'Time',A.l1_time, ...
                                             'Temp',A.ret_temp,'Z',A.ret_z(:)','Spacing',Sp);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solar zenith angle and the AIRS day/night mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SunLat,SunLon] = ndgrid(-90:10:90,-180:20:180);
SunTimes = [datenum(2000,1,1,12,0,0), datenum(2003,3,20,7,0,0), datenum(2007,6,21,12,0,0), ...
            datenum(2012,9,22,18,30,0), datenum(2016,12,21,0,0,0), datenum(2024,2,29,6,15,0)];
SunMask = NaN([size(SunLat),numel(SunTimes)]); SunSZA = SunMask;
for iT = 1:numel(SunTimes)
  Time = SunTimes(iT).*ones(size(SunLat));
  SunMask(:,:,iT) = which_airs_retrieval(SunLon,SunLat,Time);
  %the angle itself, by the same formulas as which_airs_retrieval (which returns only the mask)
  D = Time - datenum(2000,1,1,12,0,0);
  g = (357.529 + 0.98560028 .* D) .* pi ./ 180;
  q = 280.459 + 0.98564736 .* D;
  L = (q + 1.915 .* sin(g) + 0.020 .* sin(2 * g)) .* pi ./ 180;
  e = (23.439 - 0.00000036 * D) .* pi ./ 180;
  dec = asin(sin(e) .* sin(L));
  ra = atan2(cos(e) .* sin(L), cos(L));
  GMST = 18.697374558 + 24.06570982441908 .* D;
  LST = GMST + SunLon ./ 15;
  h = LST ./ 12 .* pi - ra;
  SunSZA(:,:,iT) = acos(sin(SunLat.*pi./180) .* sin(dec) + cos(SunLat.*pi./180) .* cos(dec) .* cos(h)) .* 180 / pi;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Meta.Created = datestr(now,'yyyy-mm-dd HH:MM');
Meta.Matlab  = version;
Meta.MatlabFunctionsCommit = MatlabFunctionsCommit;
Meta.Function = 'prep_airs_3d.m, which_airs_retrieval.m (DatasetSpecific/AIRS)';  %no local paths stored
save(OutFile,'Meta','DistA','DistB','DistArc','DistAz','Swath','SunLat','SunLon','SunTimes','SunMask','SunSZA','-v7');
fprintf('wrote %s\n',OutFile);
