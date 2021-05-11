%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Francesca Puppo, 2021
%
%%%%% Spike sorting
%
%
%%%%% Standard spike sorting algorithm based on PCA and k-means clustering
% Read .spk files generated from the Axion Biosystem software (must have AxIS Matlab functions to read .spk files)
% and transform electrode recordings (each electrode can record the activity of multiple neurons) into single-cell activities
%
%
% The last update in spike_sorting.m includes generation of 'raster wave' (in raster_wave, each spike waveforms
% is replaced with an unitary dirac delta). The padding.m function now has
% 3 outputs: time_pad, sp_wave and spk_wave_raster. In the connectiviy
% method, I am going to use spk_wave_raster because I observed improved
% accuracy and easier control over the parameters thanks i) less noise in
% the correlograms, ii) correlations that are now independent on the amplitude
% of the signals and iii) staigtforward definition of the minimimum correlation
% peak.
% NOTE: in the connectivity method, delta diracs corresponding to detected
% spikes are replaced with unit step functions of width 'raster_width' that
% the user can define (see main_reconstruction_MEA.m)



clear all; close all;

addpath('.\AxIS MATLAB Files')

%% Parameters

% data folder
folder = input('Input name of directory with MEA spk data to sort: ');          %set name of folder with spk data
ff = dir([folder '*spk']);

for ifile=1:size(ff,1)
    
    filename = [folder ff(ifile).name];
    
    % saving options (if 'save=1', save spike sorted data in spk_folder)
    spk_folder = [folder 'spkdata_' ff(ifile).name '\'];
    save_mode = 1;
    
    % running mode for spike sorting
    % if 0, the algorithm will autmotically spike-sort data from all wells and all electrodes;
    % if 1, the user will be able to input the specific well and electrode to check
    mode = 0;
    
    % define minimum number of spikes that an electrode must have recorded to
    % be considered an active electrode
    min_detected_spikes = 5;
    
    % plotting options
    exvariance_mode = 0;        % plot the fraction of variance explained by principal components
    cluster_plot = 0;           % visualize detected clusters in PC space
    sorted_spike_plot = 0;      % plot all spikes recorded by the eletrode under analysis, overlaid in the same graph, and visualize which one beong to which cluster
    
    
    %% Read a .spk file and load the data in Matlab
    
    FileData = AxisFile(filename);
    AllData = FileData.DataSets.LoadData;
    [num_well_rows, num_well_cols, num_eltr_rows, num_eltr_cols] = size(AllData);
    
    
    %% Choose well and electrode to analyze
    
    if mode==0
        well_R = 1:num_well_rows;
        well_C = 1:num_well_cols;
        eltr_R = 1:num_eltr_rows;
        eltr_C = 1:num_eltr_cols;
    else
        disp('Select well and electrode you want to analyze:')
        well_R = input('well row =');
        well_C = input('well colummn =');
        eltr_R = input('electrode row =');
        eltr_C = input('electrode column =');
    end
    
    
    %% Spike Sorting
    
    I_well_rows = repmat(1:length(well_R),1,length(well_C));
    J_well_colums = repmat(1:length(well_C),length(well_R),1);
    J_well_colums = J_well_colums(:);
    
    I_electrode_columns = repmat(1:length(eltr_C),1,length(eltr_R));
    J_electrode_rows = repmat(1:length(eltr_C),length(eltr_R),1);
    J_electrode_rows = J_electrode_rows(:);
    
    for well_id = 1:length(well_R)*length(well_C)
        for eltr_id = 1:length(eltr_C)*length(eltr_R)
            
            if (~isempty(AllData{well_R(I_well_rows(well_id)),well_C(J_well_colums(well_id)),eltr_C(I_electrode_columns(eltr_id)),eltr_R(J_electrode_rows(eltr_id))})) && ...
                    (size(AllData{well_R(I_well_rows(well_id)),well_C(J_well_colums(well_id)),eltr_C(I_electrode_columns(eltr_id)),eltr_R(J_electrode_rows(eltr_id))},2)>min_detected_spikes)      % discard electrodes that have recorded less than min_detected_spikes spikes
                
                [t,v] = AllData{well_R(I_well_rows(well_id)),well_C(J_well_colums(well_id)),eltr_C(I_electrode_columns(eltr_id)),eltr_R(J_electrode_rows(eltr_id))}.GetTimeVoltageVector;
                maxt = ceil(max(max(t)));       % maximum time: last time point of the latest recorded spike in a well (maximum across all electrodes in a well)
                
                % Clean data
                [preprocessed_spike_mat] = preprocess_data(v,t);
                
                
                % PCA and clustering
                [neuron_spk_V, neuron_spk_T] = cluster_data(t,v,min_detected_spikes,exvariance_mode,cluster_plot,sorted_spike_plot);
                
                
                % Padding
                emptyCells = cellfun(@isempty,neuron_spk_V);
                neuron_spk_V(emptyCells) = [];
                if size(neuron_spk_V,1)~=0
                    [spk_wave, spk_wave_raster, time_pad] = padding(neuron_spk_V, neuron_spk_T, maxt);
                    
                    
                    
                    if save_mode==1
                        if ~exist(spk_folder, 'dir')
                            mkdir(spk_folder)
                        end
                        spk_wave = cellfun(@(x) sparse(x),spk_wave,'UniformOutput',false);
                        time_pad = [time_pad(1); time_pad(end); length(time_pad)];
                        save([spk_folder 'raster_spkwave_WELL' num2str(well_R(I_well_rows(well_id))) num2str(well_C(J_well_colums(well_id)))...
                            '_EL' num2str(eltr_C(I_electrode_columns(eltr_id))) num2str(eltr_R(J_electrode_rows(eltr_id))) '.mat'],'spk_wave', 'time_pad','spk_wave_raster');
                    end
                    disp(['SUCCESS!!!! - ' num2str(size(neuron_spk_V,1)) ' neuron(s)'])
                else
                    disp('Inactive neurons');
                end
            else
                disp('No spikes')
            end
        end
        
    end
    
    
    
end
