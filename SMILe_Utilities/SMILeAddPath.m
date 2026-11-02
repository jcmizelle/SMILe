function [] = SMILeAddPath()
%%
%--------------------------------------------------------------------------
% Simple function to load all commonly used folders to the search path
%--------------------------------------------------------------------------

%%
%--------------------------------------------------------------------------
% Get a list of the current search path:
%--------------------------------------------------------------------------
fullPathList = strsplit(path(),pathsep); 

%--------------------------------------------------------------------------
%Get the path to my code directory, and exclude all (sub)directories that
% I don't want pre-loaded:
%--------------------------------------------------------------------------
basePath = SMILe_Genpath_Exclude('path/to/SMILe',...
    {'eeglab.*', 'Graph-Minimum.*', 'SMILe_Statistics_Plotting.*',...
    'SMILe_Statistics.*', 'SMILe_FT_Processing.*', 'SMILe_SourceImaging.*',...
    'SMILe_Other.*', 'SMILe_Toy.*'});

basePathList = regexp(basePath,pathsep,'Split');

basePathList = basePathList(~cellfun(@isempty,basePathList));

%--------------------------------------------------------------------------
% Loop across the list, referencing if it is already in the search path:
%--------------------------------------------------------------------------
for nfolders = 1:size(basePathList,2)
  
  if any(strcmp(fullPathList,basePathList{nfolders}))
     disp([basePathList{nfolders},' already in the path'])
  
  else 
     disp(['adding ',basePathList{nfolders},' to the path'])
     addpath(basePathList{nfolders})
  
  end
  
end

%--------------------------------------------------------------------------
% End for function:
%--------------------------------------------------------------------------
end