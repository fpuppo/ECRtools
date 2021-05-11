%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020, Francesca Puppo
%
% This function allows the user to test the connectivity method parameters 
%
% Parameters include:
% - sigma_filterig: 1 x Ns vector of standard deviations (widths) of the Gaussian filter used
% to smooth the correlation functions. It controls the generation (detection) of false
% positives (apparent and indirect connections)
% - cor_delay_interval: 1 x Nd vector of time windows over which to search for correlation
% peaks. It defines the maximum correlation delay. This parameter controls
% the filtering of false positive connections via search of correlation triangles (combinations of delays) 
% - epsilon: accepted degree of temporal approximation (a good approximation is half the width of the correlation peaks)
% NOTE: Ns and Nd can change depending on the level of granularity in the
% statistical analysis. However, epsilon could vary too. cor_delay_interval is the most critical parameter and must always vary 
% - minimum_peak: estimated hard threshold for minimum correlation peak (define the noise in the correlation function)
% - peak_discriminator_coefficient: if the peak is higher than minimun_peak, how higher than the maximum recorded peak across all
% neurons must be to be selected
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function test_parameters_updated( V, T, input_vec, sigma_test, cor_delay_test, correlation_parameters, signal_processing_parameters)

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
    if isempty(thr_spikes)||isinf(thr_spikes)||isnan(thr_spikes)||(thr_spikes<=0)
        thr_spikes = 4;
        warning('thr_spikes is set to default value 4')
    end
else
    thr_spikes = 4;
    warning('raster_width is set to default value 4')
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


%%

dt = diff(T{1}(1:2,1));

time_FFT = tic;
[aux_time1, corr1] = fast_corr(V{input_vec(1)},V{input_vec(2)}, T{input_vec(1)}, T{input_vec(2)}, dt, round(cor_delay_test/dt));
fprintf('time to evaluate FFTs = %.5g \n',toc(time_FFT))

figure(100), subplot(2,1,1),
for j =1:2
    [TT_ind, I, ~] = unique(ceil( T{input_vec(j)}(:)/dt ));
    VV = V{input_vec(j)}(:);
    VV = VV(I);
    plot(dt*(1:max(TT_ind)), sparse(TT_ind, 1 , VV )), hold on
end
N = length(T);

if effective_apply_sigma_thr==0
    % Compute correlation threshold based on sigma    
    
        % nonzero_signal_points = cellfun(@(x) size(x,2),V).';
    % sqrt_of_spikes = (min( nonzero_signal_points(:,ones(N,1)) , nonzero_signal_points(:,ones(N,1)).' ))*(floor(raster_width/2)*2+1);
    % max_peak = cellfun(@(x) max(max(abs(x),[],2)), V).';
    % minimum_peak = max_peak*max_peak.'.*sqrt_of_spikes*peak_discriminator_coeff;
    nonzero_signal_points = cellfun(@(x) size(x,2),V).';    % number of peaks in the signals
    num_spikes = min( nonzero_signal_points(:,ones(N,1)) , nonzero_signal_points(:,ones(N,1)).' );  % minimum number of peaks between two signals
    for j = input_vec(1)
        for k = input_vec(2)
            
            test_vec1 = T{j}(1,:);
            test_vec2 = T{k}(1,:);
            
            aux_mat = ones(length(test_vec2), 1) * test_vec1;
            aux_mat2 = aux_mat - test_vec2';
            
            allowed_ind = (aux_mat2 < max(cor_delay_test)) & (aux_mat2 > - max(cor_delay_test));
            N1 = sum(any(allowed_ind,1));
            N2 = sum(any(allowed_ind,2));
            
            num_spikes_time_window = min(N1,N2);
            if num_spikes_time_window>5
                if num_spikes(j,k)*.1 > num_spikes_time_window
                    minimum_peak(j,k) = inf;
                else
                    minimum_peak(j,k) = k_coeff/41*raster_width*num_spikes_time_window^alpha_coeff;
                end
            else
                minimum_peak(j,k) = inf;
                
            end
            minimum_peak(k,j) = minimum_peak(j,k);
        end
    end


else
    
    for j = input_vec(1)
        for k = input_vec(2)
            minimum_peak(j,k) = EC_sigma_thr * (std(corr1));
        end
    end
    
end


figure(100), subplot(2,1,2),hold on,plot(-aux_time1, corr1)
plot(-aux_time1, minimum_peak(input_vec(1),input_vec(2))+0*aux_time1)

for j=1:length(sigma_test)
    iF_filtered = gaussfilt(aux_time1,corr1,sigma_test(j));
    hold on, plot(-aux_time1, iF_filtered),
end




