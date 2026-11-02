function [EEG] = SMILe_REST(EEG)
%--------------------------------------------------------------------------
%
% This function re-references a dataset using the Reference Electrode 
% Standardization Technique (REST)
%
% Input:
% EEG:   Preprocesed data set from our SMILe routine. Can be epoched or 
% continuous.The original reference must be an average reference, so this 
% function call must come after the re-reference in the proprocessing 
% script.
%
%
% Output:
% EEG:   A matrix (nchan X time or nchan X time X epochs) of data 
% that have been re-referenced to "infinity."
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Make a copy of the original data:
%--------------------------------------------------------------------------
data = EEG.data;

%--------------------------------------------------------------------------
% Reformat data as needed - should be continuous for the remainder of the 
% script until declaring output:
%--------------------------------------------------------------------------
if length(size(data)) == 3
    data = reshape(data, size(data,1), size(data2)*size(data,3));
    
end

%--------------------------------------------------------------------------
% Calculate the leadfield based on a precomputed 3-concentric spheres 
% headmodel:
%--------------------------------------------------------------------------
disp('Calculating leadfield based on 3-concentric spheres headmodel.');

%--------------------------------------------------------------------------
% Use xyz coordinates from EEG.chanlocs:
%--------------------------------------------------------------------------
nchan = 1:size(data,1);
xyz_elec = zeros(length(nchan),3);
for nc = 1:length(nchan)
    xyz_elec(nc,1) = EEG.chanlocs(nchan(nc)).X;
    xyz_elec(nc,2) = EEG.chanlocs(nchan(nc)).Y;
    xyz_elec(nc,3) = EEG.chanlocs(nchan(nc)).Z;
end

%--------------------------------------------------------------------------
% Load fixed dipoles and define their oritations:
%--------------------------------------------------------------------------
[restpath, ~, ~] = fileparts(which('rest_refer.m'));
xyz_dipoles = load([restpath, filesep, 'corti869-3000dipoles.dat']);

%--------------------------------------------------------------------------
% Calculate the dipole orientations:
%--------------------------------------------------------------------------
xyz_dipOri = bsxfun (@rdivide, xyz_dipoles, sqrt(sum(xyz_dipoles.^2, 2)));
xyz_dipOri( 2601:3000, 1 ) = 0;
xyz_dipOri( 2601:3000, 2 ) = 0;
xyz_dipOri( 2601:3000, 3 ) = 1;

%--------------------------------------------------------------------------
% Define the headmodel:
%--------------------------------------------------------------------------
headmodel        = [];
headmodel.type   = 'concentricspheres';
headmodel.o      = [ 0.0000 0.0000 0.0000 ];
headmodel.r      = [ 0.8700, 0.9200, 1];
headmodel.cond   = [ 1.0000, 0.0125, 1];
headmodel.tissue = { 'brain' 'skull' 'scalp' };

%--------------------------------------------------------------------------
% Calculate and then transpose leadfield:
%--------------------------------------------------------------------------
[G, ~] = calc_leadfield3(xyz_elec, xyz_dipoles, xyz_dipOri, headmodel);

G = G';

%--------------------------------------------------------------------------
% Do the re-referencing:
%--------------------------------------------------------------------------
REST_EEG = rest_refer(data, G);

%--------------------------------------------------------------------------
% Reshape if input data were epoched:
%--------------------------------------------------------------------------
if length(size(data)) == 3
   REST_EEG = reshape(REST_EEG, size(data,1), size(data,2), size(data,3));
   
end

%--------------------------------------------------------------------------
% Insert newly re-referenced data back into the EEG structure:
%--------------------------------------------------------------------------
EEG.data = [];
EEG.data = REST_EEG;


%--------------------------------------------------------------------------
% End of function:
%--------------------------------------------------------------------------
end
%%