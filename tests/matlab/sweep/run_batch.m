% Sweep (see README.md): run one batch of get_limbsounders calls and save each
% output, or its error, as <OutDir>/<case id>.mat.
%
% Environment: SWEEP_BATCH (batch .json file), SWEEP_OUT (folder for outputs),
% ATMOSGW_MATLAB_CODE (folder of the group's MATLAB code, added with genpath).

BatchFile = getenv('SWEEP_BATCH');
OutDir    = getenv('SWEEP_OUT');
addpath(genpath(getenv('ATMOSGW_MATLAB_CODE')));

Cases = jsondecode(fileread(BatchFile));
if isstruct(Cases); Cases = num2cell(Cases); end
for iCase = 1:numel(Cases)
  C = Cases{iCase}; O = C.options;
  Args = {'Verbose',false};
  if isfield(O,'file_source') && O.file_source; Args = [Args,{'FileSource',true}];    end
  if isfield(O,'dx')          && O.dx;          Args = [Args,{'dx',true}];            end
  if isfield(O,'original_z')  && O.original_z;  Args = [Args,{'OriginalZ',true}];     end
  if isfield(O,'hindley_pws') && O.hindley_pws; Args = [Args,{'GetHindleyPWs',true}]; end
  Err = ''; D = struct(); tic;
  try
    D = get_limbsounders(C.datenums(:)', C.inst, Args{:});
  catch E
    Err = sprintf('%s: %s', E.identifier, E.message);
  end
  Seconds = toc;
  save(fullfile(OutDir,[C.id,'.mat']),'D','Err','Seconds','-v7');
end
