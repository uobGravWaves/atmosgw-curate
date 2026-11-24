%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make MATLAB reference outputs for the atmosgw.curate building blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Writes ../data/matlab_compat_fixtures.mat: small test inputs and the MATLAB
%outputs for them, compared against by tests/test_matlab_compat.py,
%tests/test_atmosphere.py and tests/test_geodesy.py. Inputs are saved with the
%outputs, so Python never has to regenerate them.
%
%Functions covered:
%  MATLAB:            interp1, unwrap, wrapToPi, deg2rad, rad2deg
%  Mapping Toolbox:   azimuth, reckon, km2deg (spherical defaults)
%  group's own code:  h2p, p2h, nph_haversine (corwin365/MatlabFunctions)
%
%How to run (set the two variables first; the run records them in the file):
%   MatlabFunctionsDir    = '<path to a checkout of corwin365/MatlabFunctions>';
%   MatlabFunctionsCommit = '<its commit hash>';
%   run('<full path>/make_compat_fixtures.m')
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
OutFile = fullfile(ThisDir,'..','data','matlab_compat_fixtures.mat');
if exist(OutFile,'file'); error('%s already exists: delete it first if you mean to remake it',OutFile); end
if ~license('test','map_toolbox'); error('Mapping Toolbox needed for azimuth, reckon and km2deg'); end
addpath(fullfile(MatlabFunctionsDir,'GeophysicalAndTime'))
rng(20261124,'twister')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interp1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%each case: sample points X, values V, query points XQ. Every method is run
%with the default extrapolation and with 'extrap'
Methods = {'linear','nearest','next','previous','pchip','spline','makima'};  %not 'cubic': cubic convolution, uniform grids only (R2024b)

Xa = cumsum(0.2 + rand(1,20));                          %ascending, uneven
Cases = struct('Name',{},'X',{},'V',{},'XQ',{});
XQa = [linspace(Xa(1)-2,Xa(end)+2,57), Xa([1 5 end])];  %inside, outside, on knots and ends
Cases(end+1) = struct('Name','ascending', 'X',Xa,         'V',randn(1,20),'XQ',XQa);
Cases(end+1) = struct('Name','descending','X',fliplr(Xa), 'V',randn(1,20),'XQ',XQa);
Cases(end+1) = struct('Name','unsorted',  'X',Xa(randperm(20)),'V',randn(1,20),'XQ',XQa);
V = randn(1,20); V(8) = NaN;
Cases(end+1) = struct('Name','nan_value', 'X',Xa,         'V',V,          'XQ',XQa);
Cases(end+1) = struct('Name','midpoints', 'X',0:4,        'V',[3 1 4 1 5],'XQ',[-0.5 0 0.5 1.5 2.5 3.5 4 4.5]);
Cases(end+1) = struct('Name','nan_query', 'X',Xa,         'V',randn(1,20),'XQ',[Xa(3) NaN Xa(4)]);
Cases(end+1) = struct('Name','three',     'X',[1 2.5 3],  'V',[2 -1 4],   'XQ',linspace(0,4,17));
%as used in get_limbsounders: a descending, uneven limb-sounder profile put onto 18:0.5:60 km
Zp = fliplr(sort(15 + 90.*rand(1,60)));
Cases(end+1) = struct('Name','profile',   'X',Zp,'V',200 + 0.8.*Zp + 5.*sin(Zp./3),'XQ',18:0.5:60);

Interp = struct();
for iC=1:numel(Cases)
  C = Cases(iC);
  S = struct('Name',C.Name,'X',C.X,'V',C.V,'XQ',C.XQ);
  for iM=1:numel(Methods)
    S.(Methods{iM})              = interp1(C.X,C.V,C.XQ,Methods{iM});
    S.([Methods{iM},'_extrap'])  = interp1(C.X,C.V,C.XQ,Methods{iM},'extrap');
  end
  Interp.(C.Name) = S;
end
clear C S iC iM V Xa XQa Zp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% unwrap, wrapToPi, deg2rad, rad2deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%unwrap: a random walk with large jumps, jumps of exactly +-pi, 2pi and 3pi and
%either side of pi, and longitudes crossing the dateline (as in get_limbsounders)
Unwrap.random   = cumsum([0, 2.5.*randn(1,200)]);
Unwrap.exact    = cumsum([0, pi, -pi, 2*pi, -2*pi, 3*pi, -3*pi, pi, pi, -pi, -pi]);
Unwrap.nearpi   = cumsum([0, pi-1e-12, pi+1e-12, -(pi-1e-12), -(pi+1e-12), 0.1]);
Unwrap.dateline = deg2rad([170 175 179.5 -179.5 -175 -170 -175 179 172]);
for F = fieldnames(Unwrap)'; UnwrapOut.(F{1}) = unwrap(Unwrap.(F{1})); end

%wrapToPi: exact multiples of pi, values just beyond, large values, NaN
WrapIn = [0, pi, -pi, 2*pi, -2*pi, 3*pi, -3*pi, 4*pi, -4*pi, pi+1e-12, -pi-1e-12, ...
          pi-1e-12, 100, -100, 1e6, 20.*randn(1,100), NaN];
WrapOut = wrapToPi(WrapIn);

%deg2rad and rad2deg
DegIn = [0, 90, 180, -180, 360, 720.*rand(1,100) - 360];
DegToRad = deg2rad(DegIn);
RadIn = [0, pi, -pi, 2*pi, 20.*randn(1,100)];
RadToDeg = rad2deg(RadIn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% h2p and p2h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%altitudes (km): negative, layer bases and either side of them, random, very high, NaN
Bases = [0 11 20 32 47 51 71];
H2pIn = [-1, -1e-9, Bases, Bases + 1e-9, Bases(2:end) - 1e-9, 130.*rand(1,100) - 5, 150, 200, NaN];
H2pOut = h2p(H2pIn);
H2pMatrix = 60.*rand(4,5); H2pMatrixOut = h2p(H2pMatrix);   %shape kept

%pressures (hPa): layer-base pressures, random log-uniform, very low, NaN,
%and a pressure above 1013.25 hPa after a valid one (p2h reuses the previous
%element's layer constants there)
BaseP = [1013.25, 226.3210, 54.7489, 8.6802, 1.1091, 0.6694, 0.0396];
P2hIn = [BaseP, BaseP.*(1+1e-9), BaseP.*(1-1e-9), 10.^(6.*rand(1,100) - 3), 1e-5, NaN, 500, 1050];
P2hOut = p2h(P2hIn);
P2hMatrix = 10.^(3.*rand(3,4)); P2hMatrixOut = p2h(P2hMatrix);
%a pressure above 1013.25 hPa as the first element: record whether MATLAB stops with an error
try;   p2h([1050, 500]); P2hFirstHighErrors = false;
catch; P2hFirstHighErrors = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% geodesy: nph_haversine, azimuth, reckon, km2deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%point pairs [lat lon]: random, identical, antipodal, across the dateline, at
%and towards the poles, due north/south/east/west
N = 100;
Pa = [180.*rand(N,1) - 90, 360.*rand(N,1) - 180];
Pb = [180.*rand(N,1) - 90, 360.*rand(N,1) - 180];
Special = [ 10  20   10   20;   %identical
            30  40  -30 -140;   %antipodal
             5 179    5 -179;   %across the dateline
            90   0   45  10;    %from the north pole
           -90   0  -45  10;    %from the south pole
            45  10   90   0;    %to the north pole
            45  10  -90   0;    %to the south pole
             0   0   10   0;    %due north
             0   0  -10   0;    %due south
             0   0    0  10;    %due east
             0   0    0 -10];   %due west
GeoA = [Pa; Special(:,1:2)];
GeoB = [Pb; Special(:,3:4)];
Haversine = nph_haversine(GeoA,GeoB);
Azimuth   = azimuth(GeoA(:,1),GeoA(:,2),GeoB(:,1),GeoB(:,2),'degrees');

%reckon (arc length and azimuth in degrees), including negative and zero arc
%lengths (HIRDLS geolocation correction uses negative ones), as a matrix as in
%module_load_HIRDLS
ReckonLat = [180.*rand(N,1) - 90;  0;  0; 89.9; -89.9; 10];
ReckonLon = [360.*rand(N,1) - 180; 179; 0;  0;    0;  -179.5];
ReckonArc = [40.*randn(N,1);        2; 0;  1;    1;    -1];
ReckonAz  = [360.*rand(N,1);       90; 45; 0;  180;   90];
[ReckonLatOut,ReckonLonOut] = reckon(ReckonLat,ReckonLon,ReckonArc,ReckonAz);
RM = struct('Lat',60.*rand(3,4) - 30,'Lon',360.*rand(3,4) - 180,'Arc',0.2.*randn(3,4),'Az',360.*rand(3,4));
[RM.LatOut,RM.LonOut] = reckon(RM.Lat,RM.Lon,RM.Arc,RM.Az);
ReckonMatrix = RM;

%km2deg (default radius, 6371 km)
KmIn = [0, 1, -1, 111.19, 20015, 1000.*randn(1,100)];
KmToDeg = km2deg(KmIn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Meta = struct('Created',datestr(now,'yyyy-mm-dd HH:MM:SS'),'Matlab',version, ...
              'MatlabFunctionsCommit',MatlabFunctionsCommit,'Seed',20261124);
save(OutFile,'Meta','Interp','Methods','Unwrap','UnwrapOut','WrapIn','WrapOut', ...
     'DegIn','DegToRad','RadIn','RadToDeg', ...
     'H2pIn','H2pOut','H2pMatrix','H2pMatrixOut','P2hIn','P2hOut','P2hMatrix','P2hMatrixOut','P2hFirstHighErrors', ...
     'GeoA','GeoB','Haversine','Azimuth', ...
     'ReckonLat','ReckonLon','ReckonArc','ReckonAz','ReckonLatOut','ReckonLonOut','ReckonMatrix', ...
     'KmIn','KmToDeg','-v7');
d = dir(OutFile); fprintf('saved %s (%.0f kB)\n',OutFile,d.bytes./1e3);
