%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020, Francesca Puppo
%
% This function reconstructs the effective connectivity of the input
% neuronal network
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [adjmat, delay_adjmat] = EC_reconstruction(V, T, sigma_filtering, cor_delay_interval, epsilon, peak_discriminator_coeff,raster_width)

N = length(V);                                        % number of neurons

nonzero_signal_points = cellfun(@(x) size(x,2),V).';    % estimate of number of peaks in the signals 
sqrt_of_spikes = sqrt(min( nonzero_signal_points(:,ones(N,1)) , nonzero_signal_points(:,ones(N,1)).' ))*(floor(raster_width/2)*2+1);  % estimate of minimum number of peaks between two signals 

minimum_peak = sqrt_of_spikes*peak_discriminator_coeff; % estimate of minimum cut-off peak

delta_t = max(cor_delay_interval)*2;                    % safe range for storing correlation vectors while saving memory and computation time

dt = diff(T{1}(1:2,1));


%% ------------------------------------------------------------------
% Compute correlations between neurons' activities 

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
            invF{n1,n2} = zeros(t_cut{n1,n2});
        end
    end
    fprintf('time to evaluate correlations for the %.5g neuron = %.5g \n',n1,toc(time_FFT1))
end
fprintf('total time to evaluate correlations = %.5g \n',toc(time_FFT))



%% ------------------------------------------------------------------
% Connectivity reconstruction

I_sigma_filtering = repmat(1:length(sigma_filtering),1,length(cor_delay_interval));
J_cor_delay_interval = repmat(1:length(cor_delay_interval),length(sigma_filtering),1);
J_cor_delay_interval = J_cor_delay_interval(:);
Z_epsilon = repmat(1:length(epsilon),1,length(J_cor_delay_interval));

J_cor_delay_interval = repmat(J_cor_delay_interval.',length(epsilon),1);
J_cor_delay_interval = J_cor_delay_interval(:);
I_sigma_filtering = repmat(I_sigma_filtering,length(epsilon),1);
I_sigma_filtering = I_sigma_filtering(:);

adjmat = cell(size(J_cor_delay_interval));
delay_adjmat = cell(size(J_cor_delay_interval));

time_adj = tic;
parfor j = 1:length(sigma_filtering)*length(cor_delay_interval)*length(epsilon)
    [adjmat{j}, delay_adjmat{j}] = adjacency_mat(epsilon(Z_epsilon(j)),cor_delay_interval(J_cor_delay_interval(j)),sigma_filtering(I_sigma_filtering(j)),N,flag_iF,invF,t_cut,minimum_peak);
end
fprintf('time to evaluate adjacency matrix = %.5g \n',toc(time_adj))



function [adjmat, delay_adjmat] = adjacency_mat(epsilon,cor_delay_interval,sigma_filtering,N,flag_iF,invF,t_cut,minimum_peak)

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
            ind = find( peaks{n1,n2}>=minimum_peak(n1,n2)); 
            peaks{n1,n2} = peaks{n1,n2}(ind);
            delays{n1,n2} = delays{n1,n2}(ind);
            delays{n2,n1} = -delays{n1,n2};
            peaks{n2,n1} = peaks{n1,n2};
        end
    end
end


% -----------------------------------------------------------------
% Detect indirect and apparent connectivity and discard correlation
% peaks based on their amplitude

[adjmat, delay_adjmat] = edge_covering_minimization(delays, peaks, epsilon);

