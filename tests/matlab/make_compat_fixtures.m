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
%  2026/11/25 [Claude] h2p/p2h inputs for the U.S. Standard Atmosphere 1976 rewrite
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

%altitudes (km): U.S. Standard Atmosphere 1976 layer bases (geometric) and either side of them, the
%lower limit and below it, 86 km (start of the table), random, the table top and above it, NaN.
%(The rand calls are unchanged from the earlier version, so the other fixtures are unchanged.)
%(layer bases and their pressures, and the Table I pressures, as computed by the 2026/11/25 h2p.m/p2h.m
%constants; written out so this script needs only h2p and p2h)
Bases = [0, 11.019067832000108, 20.06312368170136, 32.1619032229809, 47.35009222212044, 51.41247962579011, 71.80197067469582, 85.99995290624202];
BasePres = [1013.25, 226.3206397346292, 54.74888669677777, 8.680186847552283, 1.1090630555496597, 0.6693887311868731, 0.03956420428040728];
TablePres = [0.0037338, 0.0034163, 0.0031259, 0.0028602, 0.0026173, 0.0023951, ...
             0.0021919, 0.002006, 0.0018359, 0.0016804, 0.0015381, 0.0014078, ...
             0.0012887, 0.0011798, 0.0010801, 0.00098896, 0.0009056, 0.00082937, ...
             0.00075966, 0.00069592, 0.00063765, 0.00058439, 0.00053571, 0.00049122, ...
             0.00045057, 0.00041342, 0.00037948, 0.00034846, 0.00032011, 0.00027192, ...
             0.00023144, 0.00019742, 0.00016882, 0.00014477, 0.00012454, 0.00010751, ...
             9.3188e-05, 8.1142e-05, 7.1042e-05, 6.2614e-05, 5.5547e-05, 4.957e-05, ...
             4.4473e-05, 4.0096e-05, 3.6312e-05, 3.3022e-05, 3.0144e-05, 2.7615e-05, ...
             2.5382e-05, 2.3401e-05, 2.1635e-05, 2.0055e-05, 1.8635e-05, 1.7354e-05, ...
             1.6194e-05, 1.5141e-05, 1.4183e-05, 1.3307e-05, 1.2505e-05, 1.1769e-05, ...
             1.1092e-05, 1.0468e-05, 9.8907e-06, 9.3568e-06, 8.8617e-06, 8.4018e-06, ...
             7.9739e-06, 7.5751e-06, 7.2028e-06, 6.855e-06, 6.5294e-06, 6.2243e-06, ...
             5.938e-06, 5.6691e-06, 5.4162e-06, 5.1781e-06, 4.9538e-06, 4.7421e-06];
H2pIn = [-5, -5-1e-9, -1, -1e-9, Bases, Bases + 1e-9, Bases(2:end) - 1e-9, 130.*rand(1,100) - 5, 86, 100, 149, 150, 200, NaN];
H2pOut = h2p(H2pIn);
H2pMatrix = 60.*rand(4,5); H2pMatrixOut = h2p(H2pMatrix);   %shape kept

%pressures (hPa): layer-base pressures and either side of them, Table I pressures, random log-uniform,
%the limits and beyond them, above 1013.25 hPa, zero, negative, NaN
P2hIn = [BasePres, BasePres.*(1+1e-9), BasePres.*(1-1e-9), 10.^(6.*rand(1,100) - 3), TablePres, 1777.6, 1800, 1050, 4.7e-6, 0, -1, NaN];
P2hOut = p2h(P2hIn);
P2hMatrix = 10.^(3.*rand(3,4)); P2hMatrixOut = p2h(P2hMatrix);

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
     'H2pIn','H2pOut','H2pMatrix','H2pMatrixOut','P2hIn','P2hOut','P2hMatrix','P2hMatrixOut', ...
     'GeoA','GeoB','Haversine','Azimuth', ...
     'ReckonLat','ReckonLon','ReckonArc','ReckonAz','ReckonLatOut','ReckonLonOut','ReckonMatrix', ...
     'KmIn','KmToDeg','-v7');
d = dir(OutFile); fprintf('saved %s (%.0f kB)\n',OutFile,d.bytes./1e3);
