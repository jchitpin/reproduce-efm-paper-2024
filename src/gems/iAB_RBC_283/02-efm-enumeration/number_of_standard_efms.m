%% DESCRIPTION
% This script does the following:
% (1) Computes the number of (standard) EFMs in the iAB_RBC_283 network

%% USER PARAMETERS
% Set working directory of this script. Example:
cd('/home/jchitpin/Documents/PhD/Projects/reproduce-efm-paper-2024/src/gems/iAB_RBC_283/')
filename = '../../../data/gems/iAB_RBC_283/processed/stoichiometry-matrix-processed.csv';

% Path to FluxModeCalculator
addpath(genpath('/home/jchitpin/Documents/PhD/Code/MATLAB'))
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
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
S = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;

% Calculate binary EFMs (there are 5,058,711 EFMs)
rev = logical(zeros(1, size(S,2)));
[efm_bin,S_unc,id_unc,T,stats] = calculate_flux_modes(S, rev, 'MaxThreads', 1);  % 1303.4795 s
[efm_bin,S_unc,id_unc,T,stats] = calculate_flux_modes(S, rev, 'MaxThreads', 2);  %  711.7735 s
[efm_bin,S_unc,id_unc,T,stats] = calculate_flux_modes(S, rev, 'MaxThreads', 4);  %  386.9873 s
[efm_bin,S_unc,id_unc,T,stats] = calculate_flux_modes(S, rev, 'MaxThreads', 8);  %  233.2888 s
[efm_bin,S_unc,id_unc,T,stats] = calculate_flux_modes(S, rev, 'MaxThreads', 16); %  153.2752 s
[efm_bin,S_unc,id_unc,T,stats] = calculate_flux_modes(S, rev, 'MaxThreads', 32); %  161.1707 s

% Calculate coefficients of binary EFMs and check for consistency
%[efm,err]=bin2num_flux_modes(efm_bin,S);
