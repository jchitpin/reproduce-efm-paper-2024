%% DESCRIPTION
% This script does the following:
% (1) Computes the number of (standard) EFMs in the e_coli_core network.
% NOTE: FluxModeCalculator must be installed in MATLAB for this script to run!
% You may want to run the following code after manually importing the
% stoichiometry matrix.

%% USER PARAMETERS
% Set working directory of this script. Example:
cd('/home/jchitpin/Documents/PhD/Projects/reproduce-efm-paper-2024/src/gems/e_coli_core/')
filename = '../../../data/gems/e_coli_core/processed/stoichiometry-matrix-processed.csv';

% Path to FluxModeCalculator
addpath(genpath('/home/jchitpin/Documents/PhD/Code/MATLAB')) % your path to FluxModeCalculator
%% FluxModeCalculator - Ensure this MATLAB package is installed and test it with the following example
clc

% Generate random stoichiometric matrix
n=[10 20];
S=[diag(floor(rand(n(1),1)*3)+1),round(0.7*randn(n(1),n(2)-n(1)))];
A=1;
while rank(A)<n(1)
    A=round(randn(n(1)*1.5,n(1))*0.3);
end
S=A*S;
rev=rand(1,n(2))>0.6;

% Show help of both functions
help calculate_flux_modes
help bin2num_flux_modes

% Calculate binary EFMs
[efm_bin,S_unc,id_unc,T,stats]=calculate_flux_modes(S,rev);

% Calculate coefficients of binary EFMs and check for consistency
[efm,err]=bin2num_flux_modes(efm_bin,S);

%% Import final stoichiometry matrix
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
S = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;

% Calculate binary EFMs (there are 5,999,302 EFMs)
rev = logical(zeros(1, size(S,2)));
[efm_bin,S_unc,id_unc,T,stats] = calculate_flux_modes(S, rev, 'MaxThreads', 1);  % 5859.5641 s
[efm_bin,S_unc,id_unc,T,stats] = calculate_flux_modes(S, rev, 'MaxThreads', 2);  % 3337.9084 s
[efm_bin,S_unc,id_unc,T,stats] = calculate_flux_modes(S, rev, 'MaxThreads', 4);  % 1810.2888 s
[efm_bin,S_unc,id_unc,T,stats] = calculate_flux_modes(S, rev, 'MaxThreads', 8);  % 1092.3322 s
[efm_bin,S_unc,id_unc,T,stats] = calculate_flux_modes(S, rev, 'MaxThreads', 16); %  717.5590 s
[efm_bin,S_unc,id_unc,T,stats] = calculate_flux_modes(S, rev, 'MaxThreads', 32); %  772.3231 s

% Calculate coefficients of binary EFMs and check for consistency
%[efm,err]=bin2num_flux_modes(efm_bin,S);
