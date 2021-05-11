%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021, Francesca Puppo
%
% The user can:
% - define the experiment folder
% - test the connectivity parameters by using the test_parameters.m function
% - call function to load, check and correct spike sorted data (remapping_time.m)
% - run the connectivity algorithm
%
% INPUT DATA: sorry! this release still works only on MEA data formatted in
% a specific way. Next release will be more general. Please see description of data format here below.
% Contact me at fpuppo@ucsd.edu for further input
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
% - EC_delay_adjmat: inferred NxN matrix (N total number of neurons) of time
% delays between neurons connected via direct and causal links
% - frequency_matrix: NxN ajacency matrix with connection frequencies
% - EC_adjmat: inferred NxN ajacency matrix indicating effective (direct and causal) connectivity between neurons
%
% Parameters:
% 1) Signal processing parameters (preparation of input data)
% - raster_width: each spike is replaced by an unitary step of width (number of points) 'raster_width'
% - thr_spikes: minimum number of spikes of neuron j that must fall close to each spikes of neuron k to consider an existing connection  
% 2) Main parameters for reconstruction of the EC connectivity  matrix:
% - sigma_filtering: 1 x Ns vector of standard deviations (widths) of the Gaussian filter used
% to smooth the correlation functions.
% - cor_delay_interval: 1 x Nd vector of time windows over which to search for correlation
% peaks. It defines the maximum correlation delay.
% - epsilon: accepted degree of temporal approximation
% NOTE: Ns and Nd can change depending on the level of granularity in the
% statistical analysis. However, epsilon could vary too. cor_delay_interval is the most critical parameter and must always vary
% - pearson: use cross-ccorelation (when 0), or Pearson correlation (when 1)
% - discrimination_threshold
% 3) Correlation parameters (identification of correlation peaks):
% - k_coeff and alpha_coeff define the noise threshold in the correlations
% 4) If parallel computing, set the number of workers (NumWorkers) you want to use
%
%
%
% Reference reading :
% F. Puppo, D. Pre', A. Bang, G. Silva
% "Super-selective reconstruction of causal and direct connectivity with
% application to in-vitro iPSC neuronal networks"
% bioRxiv, https://doi.org/10.1101/2020.04.28.067124
%
% For support, please contact Francesca Puppo (fpuppo@ucsd.edu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all


%% =======================================================================
%%%%%%%%%%%%%%%%%%Data folder and file name%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder = ['.\spkdata\'];  % set folder

well = [21];  %set well/wells to analyze

%% =======================================================================
%%%%%%%%%%%%%%%%%%reconstruction parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parallel pool parameters
NumWorkers = 8;

% signal processing parameters
signal_processing_parameters.raster_width = 5;                             % Each spike is replaced by an unitary step of width (number of points) 'raster_width'
signal_processing_parameters.thr_spikes = 4;                               % Minimum number of spikes of neuron j that must fall close to each spikes of neuron k to consider an existing connection
signal_processing_parameters.spike_percentage = 0;

% effective connectivity reconstruction parameters
EC_reconstruction_parameters.sigma_filtering = linspace(.0001,.0005,5); %linspace(.000075,.00025,10);%[.0003 .00075 .0012];
EC_reconstruction_parameters.cor_delay_interval = linspace(10e-3,30e-3,4);%[15e-3 20e-3 25e-3];
EC_reconstruction_parameters.epsilon = 1e-3; %linspace(0.5e-3,2e-3,10);
EC_reconstruction_parameters.pearson = 1;                               % use Pearson coefficient or regular cross-correlation in triangle analysis
discrimination_threshold = 0.85;

% correlation parameters
correlation_parameters.functional = 0;                               %1: extract both functional and effective connectivity; 0: only extract effective
correlation_parameters.functional_apply_sigma_thr = 0;               %functional connetome - 0: apply standard k_coeff and alpha_coeff threshold; 1: apply sigma threshold
correlation_parameters.effective_apply_sigma_thr = 0;                %effective connectome - 0: apply standard k_coeff and alpha_coeff threshold; 1: apply sigma threshold
correlation_parameters.EC_sigma_thr = 0;                             %define sigma for sigma-thresholding in EC reconstruction
correlation_parameters.FC_sigma_thr = 0;                             %define sigma for sigma-thresholding in FC reconstruction; if 0, no thresholding will be applied (all peaks selected)
correlation_parameters.k_coeff = 7;%linspace(0,40,50);
correlation_parameters.alpha_coeff = .85;


%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%test parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uncomment to test specific parameters

% % parameters to test
% test_neurons =[3 25];
% test_sigma = EC_reconstruction_parameters.sigma_filtering([1,end]);
% test_cor_delay = 30e-3;
% h = 31;                                 % set well to analyze
% folder = [baseDir folders(1).name '\'];
% if ~isempty(dir([folder '*raster_spkwave_WELL*' num2str(h) '_*']))
%
%     % load MEA spike sorted data and remap time
%     [time, V, T] = remapping_time(folder,well,signal_processing_parameters,EC_reconstruction_parameters.cor_delay_interval, MEA);
%
%     % test parameters
%     test_parameters_updated(V, T, test_neurons,test_sigma , test_cor_delay, correlation_parameters, signal_processing_parameters)
% else
%     disp('no activity in this well!')
% end



%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%CONNECTIVITY ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(['*******Start analysis********'])

files = dir([folder '*raster_spkwave_WELL*' num2str(well) '_*']);
if ~isempty(files)
    
    %% load MEA spike sorted data and remap time
    [time, V, T] = remapping_time(folder,well,signal_processing_parameters,EC_reconstruction_parameters);
    
    %% paralelizzation
    if isempty(gcp('nocreate'))
        parpool(NumWorkers)
    else
        pool = gcp;
        if pool.NumWorkers~= NumWorkers
            delete(gcp('nocreate'))
            parpool(NumWorkers)
        end
    end
    
    %% run effective/functional connectivity analysis
    disp('........................................................')
    disp('..................ANALYSIS STARTED......................')
    disp('........................................................')
    [adjmat, delay_adjmat, FC_adjmat, FC_delays, FC_peaks, I_sigma_filtering, J_cor_delay_interval, Z_epsilon, W_k_coeff] = connectivity_reconstruction(V, T, EC_reconstruction_parameters, correlation_parameters, signal_processing_parameters, time);
    
    %% reconstruct adjacency matrix
    % looking for connections that always exist (for each varying sigma_filtering and cor_delay_interval values)
    frequency_matrix = 0;
    EC_delay_adjmat = inf;
    for j=1:length(adjmat)
        frequency_matrix = frequency_matrix+adjmat{j};
        EC_delay_adjmat = min(delay_adjmat{j},EC_delay_adjmat);
    end
    EC_delay_adjmat(frequency_matrix<discrimination_threshold*length(adjmat)) = inf;
    frequency_matrix(frequency_matrix<discrimination_threshold*length(adjmat)) = 0;
    EC_adjmat = frequency_matrix;
    EC_adjmat(EC_adjmat~=0) = 1;
    
    
end












