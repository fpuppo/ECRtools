%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020, Francesca Puppo
%
% Super-selective reconstruction of effective connectivity in neuronal networks
%
% The user can:
% - define the experiment folder
% - set connectivity parameters by using the test_parameters.m function
% - call function to load, check and correct spike sorted data (remapping_time.m)
% - run the connectivity algorithm (EC_reconstruction.m)
%
% Parameters:
% 1) Parameters for reconstruction of the connectivity  matrix:
% - sigma_filterig: 1 x Ns vector of standard deviations (widths) of the Gaussian filter used
% to smooth the correlation functions.
% - cor_delay_interval: 1 x Nd vector of time windows over which to search for correlation
% peaks. It defines the maximum correlation delay.
% - epsilon: accepted degree of temporal approximation
% NOTE: Ns and Nd can change depending on the level of granularity in the
% statistical analysis. However, epsilon could vary too. cor_delay_interval is the most critical parameter and must always vary
% - peak_discriminator_coefficient: it defines the minimum correlation peak (used in the spike search algorithm).
% 2) If parallel computing, set the number of workers (NumWorkers) you want to use
%
%
% The variable folder contains the path to the spike sorted data (see examples in shared folder \spkdata\).
% The spkdata folder includes N_eltr' x N_well' mat files where N_eltr' is the number of electrodes
% that have recorded some activity and N_well' is the number of wells with
% some activity. Each file in the folder is named after the MEA well and the electrode ('spkwave_WELL22_EL12.mat')
% Each file must contain the following variables:
%   -- time_pad: a 3 elements vector such that linspace(time_pad(1), time_pad(2), time_pad(3)) defines the recorded time series (time_pad(1) must be >=0)
%   -- raster_spk_wave: 1 x N cell array with N number of neurons detected by an electrode;
%   each cell is a time_pad(3) x 1 array containing all spikes of a neuron concatenated via zero-padding 
%   (raster plots = 1 in correspondence of the spike maximum; 0 everywhere else)
%
%
% Output:
% - delay_adjmat_: inferred NxN matrix (N total number of neurons) of time
% delays between neurons connected via direct and causal links
% - adjmat_: inferred NxN ajacency matrix indicating effective (direct and causal) connectivity between neurons
%
% Reference reading :
% F. Puppo, D. Pre, A. Bang, G. Silva
% "Super-selective reconstruction of causal and direct connectivity with
% application to in-vitro iPSC neuronal networks"
% bioRxiv, https://doi.org/10.1101/2020.04.28.067124
%
% For support, please contact Francesca Puppo (fpuppo@ucsd.edu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% -----------------------------------------------------------------
% set parameters

% parallel pool parameters
NumWorkers = 16;

% signal processing parameters
padding_tollerance = 1e-6;
raster_width = 41; % Each spike is replaced by an unitary step of width (number of points) 'raster_width'

% reconstruction parameters
sigma_filtering = [.0003 .00055 .0008];
cor_delay_interval = [20e-3 25e-3 30e-3];
epsilon  = 1e-3;
peak_discriminator_coeff = .25;

%% -----------------------------------------------------------------
% % uncomment to test specific parameters
%
% % parameters to test
% test_neurons =[3 10];
% test_sigma = [.0003 .00055 .0008];
% test_cor_delay = 30e-3;
% h = 21;       % set well to analyze
% folder = ['./spkdata/'];
% if ~isempty(dir([folder '*raster_spkwave_WELL*' num2str(h) '_*']))
%     
%     % load MEA spike sorted data and remap time
%     [time, V, T] = remapping_time(folder,h,raster_width,cor_delay_interval);
%     
%     % test parameters
%     test_parameters(V, T, test_neurons,test_sigma , test_cor_delay, peak_discriminator_coeff,raster_width)
% end


%% -----------------------------------------------------------------
% start analysis

h = 21;       % set well to analyze
disp(['WELL ' num2str(h)])

folder = ['./spkdata/'];  % set folder

if ~isempty(dir([folder '*raster_spkwave_WELL*' num2str(h) '_*']))
    
    % load MEA spike sorted data and remap time
    [time, V, T] = remapping_time(folder,h,raster_width,cor_delay_interval);
    
    % paralelizzation
    if isempty(gcp('nocreate'))
        parpool(NumWorkers)
    else
        pool = gcp;
        if pool.NumWorkers~= NumWorkers
            delete(gcp('nocreate'))
            parpool(NumWorkers)
        end
    end
    
    % run effective connectivity analysis
    [adjmat, delay_adjmat] = EC_reconstruction(V, T,  sigma_filtering, cor_delay_interval, epsilon, peak_discriminator_coeff, raster_width);
    
        % reconstruct adjacency matrix looking for connection that always exist
    % (for each varying sigma_filtering and cor_delay_interval values)
    adjmat_ = 0;
    delay_adjmat_ = inf;
    for j=1:length(adjmat)
        adjmat_ = adjmat_+adjmat{j};
        delay_adjmat_ = min(delay_adjmat{j},delay_adjmat_);
    end
    delay_adjmat_(adjmat_<.9*length(adjmat)) = inf;
    adjmat_(adjmat_<.9*length(adjmat)) = 0;
    
end
