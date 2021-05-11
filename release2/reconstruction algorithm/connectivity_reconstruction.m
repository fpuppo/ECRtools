%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021, Francesca Puppo
%
% This function reconstructs both the functional and effective connectivity of the input
% neuronal network
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [adjmat, delay_adjmat, FC_adjmat, FC_delays, FC_peaks, I_sigma_filtering, J_cor_delay_interval, Z_epsilon, W_k_coeff] = connectivity_reconstruction(V, T, EC_reconstruction_parameters, correlation_parameters, signal_processing_parameters, time)


if isfield(signal_processing_parameters,'raster_width')
    raster_width = signal_processing_parameters.raster_width;
    if isempty(raster_width)||isinf(raster_width)||isnan(raster_width)||(raster_width<=0)
        raster_width = 5;
        warning('raster_width is set to default value 5')
    end
    if (ceil(raster_width)~=raster_width)||(mod(ceil(raster_width),2)~=1)
        raster_width = ceil(raster_width)+mod(ceil(raster_width),2);
        warning('raster_width is set to the next odd value')
    end
else
    raster_width = 5;
    warning('raster_width is set to default value 5')
end

if isfield(signal_processing_parameters,'thr_spikes')
    thr_spikes = signal_processing_parameters.thr_spikes;
    if isempty(thr_spikes)||isinf(thr_spikes)||isnan(thr_spikes)||(thr_spikes<0)
        thr_spikes = 4;
        warning('thr_spikes is set to default value 4')
    end
else
    thr_spikes = 4;
    warning('raster_width is set to default value 4')
end


if isfield(signal_processing_parameters,'spike_percentage')
    spike_percentage = signal_processing_parameters.spike_percentage;
    if isempty(spike_percentage)||isinf(spike_percentage)||isnan(spike_percentage)||(spike_percentage<0)
        spike_percentage = 0;
        warning('spike_percentage is set to default value 0')
    end
else
    spike_percentage = 0;
    warning('spike_percentage is set to default value 0')
end


if isfield(EC_reconstruction_parameters,'sigma_filtering')
    sigma_filtering = EC_reconstruction_parameters.sigma_filtering;
    if isempty(sigma_filtering)
        sigma_filtering = linspace(.0001,.0005,5);
        warning('sigma_filtering set to its default value')
    end
else
    sigma_filtering = linspace(.0001,.0005,5);
    warning('sigma_filtering set to its default value')
end

if isfield(EC_reconstruction_parameters,'cor_delay_interval')
    cor_delay_interval = EC_reconstruction_parameters.cor_delay_interval;
    if isempty(cor_delay_interval)
        cor_delay_interval = linspace(10e-3,30e-3,4);
        warning('cor_delay_interval set to its default value')
    end
else
    cor_delay_interval = linspace(10e-3,30e-3,4);
        warning('cor_delay_interval set to its default value')
end


if isfield(EC_reconstruction_parameters,'epsilon')
    epsilon = EC_reconstruction_parameters.epsilon;
    if isempty(epsilon)
        epsilon = 1e-3;
        warning('epsilon set to its default value')
    end
else
    epsilon = 1e-3;
    warning('epsilon set to its default value')
end

if isfield(EC_reconstruction_parameters,'pearson')
    pearson = EC_reconstruction_parameters.pearson;
    if isempty(pearson)
        pearson = 1;
        warning('pearson set to its default value')
    end
else
    pearson = 1;
    warning('pearson set to its default value')
end

if isfield(correlation_parameters,'functional')
    functional = correlation_parameters.functional;
    if isempty(functional)
        functional = 1;
        warning('functional set to its default value')
    end
else
    functional = 1;
    warning('functional set to its default value')
end


if isfield(correlation_parameters,'functional_apply_sigma_thr')
    functional_apply_sigma_thr = correlation_parameters.functional_apply_sigma_thr;
    if isempty(functional_apply_sigma_thr)
        functional_apply_sigma_thr = 0;
        warning('functional_apply_sigma_thr set to its default value')
    end
else
    functional_apply_sigma_thr = 0;
    warning('functional_apply_sigma_thr set to its default value')
end

if isfield(correlation_parameters,'effective_apply_sigma_thr')
    effective_apply_sigma_thr = correlation_parameters.effective_apply_sigma_thr;
    if isempty(effective_apply_sigma_thr)
        effective_apply_sigma_thr = 0;
        warning('effective_apply_sigma_thr set to its default value')
    end
else
    effective_apply_sigma_thr = 0;
    warning('effective_apply_sigma_thr set to its default value')
end

if isfield(correlation_parameters,'EC_sigma_thr')
    EC_sigma_thr = correlation_parameters.EC_sigma_thr;
    if isempty(EC_sigma_thr)
        EC_sigma_thr = 2;
        warning('EC_sigma_thr set to its default value')
    end
else
    EC_sigma_thr = 2;
    warning('EC_sigma_thr set to its default value')
end


if isfield(correlation_parameters,'FC_sigma_thr')
    FC_sigma_thr = correlation_parameters.FC_sigma_thr;
    if isempty(FC_sigma_thr)
        FC_sigma_thr = 2;
        warning('FC_sigma_thr set to its default value')
    end
else
    FC_sigma_thr = 2;
    warning('FC_sigma_thr set to its default value')
end


if isfield(correlation_parameters,'k_coeff')
    k_coeff = correlation_parameters.k_coeff;
    if isempty(k_coeff)
        k_coeff = 7;
        warning('k_coeff set to its default value')
    end
else
    k_coeff = 7;
    warning('k_coeff set to its default value')
end

if isfield(correlation_parameters,'alpha_coeff')
    alpha_coeff = correlation_parameters.alpha_coeff;
    if isempty(alpha_coeff)
        alpha_coeff = 0.85;
        warning('alpha_coeff set to its default value')
    end
else
    alpha_coeff = 0.85;
    warning('alpha_coeff set to its default value')
end

%%


N = length(V);                                        % number of neurons


%% ------------------------------------------------------------------
% Compute correlations between neurons' activities
delta_t = max(cor_delay_interval)*2;                    % safe range for storing correlation vectors while saving memory and computation time

dt = diff(T{1}(1:2,1));

disp('************Computation of correlations**************')
time_FFT = tic;
for n1 = 1:N
    time_FFT1 = tic;
    appV = V{n1};
    appT = T{n1};
    parfor n2 = n1+1:N
        if (~isempty(appV))&&(~isempty(V{n2}))
            flag_iF(n1,n2) = true;
            if size(appV,2)>size(V{n2},2)
                [t_cut{n1,n2}, invF{n1,n2}] = fast_corr(appV,V{n2}, appT, T{n2}, dt, round(delta_t/dt));
            else
                [aux_time2, corr2] = fast_corr(V{n2},appV, T{n2}, appT,  dt, round(delta_t/dt));
                invF{n1,n2} = corr2(end:-1:1);
                t_cut{n1,n2} = -aux_time2(end:-1:1);
            end
        else
            flag_iF(n1,n2) = false;
            t_cut{n1,n2} = (-round(delta_t/dt):round(delta_t/dt)).'*dt;
            invF{n1,n2} = zeros(size(t_cut{n1,n2}));
        end
    end
    fprintf('time to evaluate correlations for the %.5g neuron = %.5g \n',n1,toc(time_FFT1))
end
fprintf('total time to evaluate correlations = %.5g \n',toc(time_FFT))



%% ------------------------------------------------------------------
% Compute mean and standard deviation of all signals

mu = zeros(N,1);
std_dev = zeros(N,1);
for n = 1:N
    raster = sparse(round(T{n}/dt),1,1,length(time),1);
    mu(n) = mean(raster);
    std_dev(n) = std(raster);
end



%% ------------------------------------------------------------------
% Compute correlation threshold and reconstruct the effective and
% functional connectome (functional = 1) or only the effective connectome
% (functional = 0)


% Compute correlation threshold for effective connectome
% Compute correlation threshold based on sigma
% Compute correlation threshold based on number of spikes of one neuron within the time window of the other neuron
[EC_minimum_peak] = compute_generalized_threshold(invF, V, T, N, EC_sigma_thr, k_coeff, alpha_coeff, cor_delay_interval, thr_spikes, raster_width, spike_percentage,effective_apply_sigma_thr);


% Reconstruct effective connectivity
disp('************Inference of effective connectome**************')

I_sigma_filtering = repmat(1:length(sigma_filtering),1,length(cor_delay_interval));
J_cor_delay_interval = repmat(1:length(cor_delay_interval),length(sigma_filtering),1);
J_cor_delay_interval = J_cor_delay_interval(:);
Z_epsilon = repmat(1:length(epsilon),1,length(J_cor_delay_interval));

J_cor_delay_interval = repmat(J_cor_delay_interval.',length(epsilon),1);
J_cor_delay_interval = J_cor_delay_interval(:);
I_sigma_filtering = repmat(I_sigma_filtering,length(epsilon),1);
I_sigma_filtering = I_sigma_filtering(:);
W_k_coeff = repmat(1:length(k_coeff),1,length(J_cor_delay_interval));

Z_epsilon = repmat(Z_epsilon,1,length(k_coeff));
Z_epsilon = Z_epsilon(:);
J_cor_delay_interval = repmat(J_cor_delay_interval.',length(k_coeff),1);
J_cor_delay_interval = J_cor_delay_interval(:);
I_sigma_filtering = repmat(I_sigma_filtering.',length(k_coeff),1);
I_sigma_filtering = I_sigma_filtering(:);

adjmat = cell(size(J_cor_delay_interval));
delay_adjmat = cell(size(J_cor_delay_interval));

time_adj = tic;
parfor j = 1:length(sigma_filtering)*length(cor_delay_interval)*length(epsilon)*length(k_coeff)
    [adjmat{j}, delay_adjmat{j}] = adjacency_mat(epsilon(Z_epsilon(j)),cor_delay_interval(J_cor_delay_interval(j)),sigma_filtering(I_sigma_filtering(j)),N,flag_iF,invF,t_cut,squeeze(EC_minimum_peak(:,:,W_k_coeff(j))),mu,std_dev,pearson);
end
fprintf('time to evaluate adjacency matrix = %.5g \n',toc(time_adj))


if functional
       
    % Compute correlation threshold for functional connectome
    % Compute correlation threshold based on sigma for functional connectome
    % Compute correlation threshold based on number of spikes of one neuron within the time window of the other neuron
    [FC_minimum_peak] = compute_generalized_threshold(invF, V, T, N, FC_sigma_thr, k_coeff, alpha_coeff, cor_delay_interval, thr_spikes, raster_width, spike_percentage,effective_apply_sigma_thr);
    
    % Detect functional links (functional connectome)
    disp('************Inference of functional connectome**************')
    [FC_delays,FC_peaks] = process_peak(N,flag_iF,invF,t_cut,min(sigma_filtering),max(cor_delay_interval),FC_minimum_peak,mu,std_dev,0);
    FC_peaks = cellfun(@(x,y) x(y<0), FC_peaks, FC_delays, 'UniformOutput', false);
    FC_delays = cellfun(@(x,y) x(y<0), FC_delays, FC_delays, 'UniformOutput', false);
    FC_adjmat = cellfun(@(x) ~isempty(x), FC_peaks);
    
else
    
    FC_peaks = [];
    FC_delays = [];
    FC_adjmat = [];
    
end



%% -------------------------------------------------------------------
% FUNCTIONS


function [minimum_peak] = compute_generalized_threshold(invF, V, T, N, sigma_thr, k_coeff, alpha_coeff, cor_delay_interval, thr_spikes, raster_width, spike_percentage, effective_apply_sigma_thr)

nonzero_signal_points = cellfun(@(x) size(x,2),V).';                                           % number of peaks in the signals
num_spikes = min( nonzero_signal_points(:,ones(N,1)) , nonzero_signal_points(:,ones(N,1)).' );  % minimum number of peaks between two signals

minimum_peak = zeros(N,N,length(k_coeff));
num_spikes_time_window = zeros(N,N);
for j = 1:N-1
    for k = j+1:N
        
        test_vec1 = T{j}(1,:);
        test_vec2 = T{k}(1,:);
        
        aux_mat = ones(length(test_vec2), 1) * test_vec1;
        aux_mat2 = aux_mat - test_vec2';
        
        allowed_ind = (aux_mat2 < max(cor_delay_interval)) & (aux_mat2 > - max(cor_delay_interval));
        N1 = sum(any(allowed_ind,1));
        N2 = sum(any(allowed_ind,2));
        
        num_spikes_time_window(j,k) = min(N1,N2);
        for h=1:length(k_coeff)
            if num_spikes_time_window(j,k)>thr_spikes
                if num_spikes(j,k)*spike_percentage > num_spikes_time_window(j,k)
                    minimum_peak(j,k,h) = inf;
                else
                    if effective_apply_sigma_thr 
                        minimum_peak(j,k) = sigma_thr * (std(invF{j,k}));
                    else
                        minimum_peak(j,k,h) = k_coeff(h)/41*raster_width*num_spikes_time_window(j,k)^alpha_coeff;
                    end
                end
                minimum_peak(k,j,h) = minimum_peak(j,k,h);
            else
                minimum_peak(j,k,h) = inf;
            end
        end
    end
end




function [adjmat, delay_adjmat] = adjacency_mat(epsilon,cor_delay_interval,sigma_filtering,N,flag_iF,invF,t_cut,minimum_peak,mu,std_dev,pearson)

% -----------------------------------------------------------------
% Apply Gaussian filtering to smooth correlation signals and find correlation peaks/delays

[delays,peaks] = process_peak(N,flag_iF,invF,t_cut,sigma_filtering,cor_delay_interval,minimum_peak,mu,std_dev,pearson);


% -----------------------------------------------------------------
% Detect indirect and apparent connectivity and discard correlation
% peaks based on their amplitude

[adjmat, delay_adjmat] = edge_covering_minimization(delays, peaks, epsilon);





function [delays,peaks] = process_peak(N,flag_iF,invF,t_cut,sigma_filtering,cor_delay_interval,minimum_peak,mu,std_dev,pearson)

% -----------------------------------------------------------------
% Apply Gaussian filtering to smooth correlation signals and find correlation peaks/delays

delays = cell(N,N);
peaks = cell(N,N);
for n1 = 1:N
    for n2 = n1+1:N
        if flag_iF(n1,n2)
            % peak search
            iF = invF{n1,n2};
            iF_filtered = gaussfilt(t_cut{n1,n2},iF,sigma_filtering);
            [~, locs] = findpeaks(iF_filtered);
            locs = locs((t_cut{n1,n2}(locs)<cor_delay_interval)&(t_cut{n1,n2}(locs)>-cor_delay_interval));
            
            % peak location refinement within 3 points distance (this can become an advanced setting)
            locs0 = locs;
            for j=-3:3
                locs1 = max(min(locs0+j,length(iF)),1);
                idx = find(iF(locs1)>iF(locs));
                locs(idx) = locs1(idx);
            end
            delays{n1,n2} = t_cut{n1,n2}(locs);
            peaks{n1,n2} = iF(locs);
            delays{n2,n1} = -t_cut{n1,n2}(locs);
            peaks{n2,n1} = iF(locs);
        else
            delays{n1,n2} = zeros(1,0);
            peaks{n1,n2} = zeros(1,0);
            delays{n2,n1} = zeros(1,0);
            peaks{n2,n1} = zeros(1,0);
        end
    end
end

% eliminate peaks smaller than the minimum
for n1 = 1:N
    for n2 = n1+1:N
        if flag_iF(n1,n2)
            % extract correlation peaks for effective connectivity
            ind = find( peaks{n1,n2}>=minimum_peak(n1,n2));
            if pearson
                peaks{n1,n2} = (peaks{n1,n2}(ind)-(mu(n1)*mu(n2)))/(std_dev(n1)*std_dev(n2));
                delays{n1,n2} = delays{n1,n2}(ind);
                delays{n2,n1} = -delays{n1,n2};
                peaks{n2,n1} = peaks{n1,n2};
            else
                peaks{n1,n2} = peaks{n1,n2}(ind);
                delays{n1,n2} = delays{n1,n2}(ind);
                delays{n2,n1} = -delays{n1,n2};
                peaks{n2,n1} = peaks{n1,n2};
            end
        end
    end
end

