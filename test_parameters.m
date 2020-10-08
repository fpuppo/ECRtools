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


function test_parameters( V, T, input_vec, sigma_test, cor_delay_test,peak_discriminator_coeff,raster_width)

dt = diff(T{1}(1:2,1));

figure(100), subplot(2,1,1),
for j =1:2
[TT_ind, I, ~] = unique(ceil( T{input_vec(j)}(:)/dt ));
VV = V{input_vec(j)}(:);
VV = VV(I);
plot(dt*(1:max(TT_ind)), sparse(TT_ind, 1 , VV )), hold on
end
N = length(T);
nonzero_signal_points = cellfun(@(x) size(x,2),V).';    
sqrt_of_spikes = (min( nonzero_signal_points(:,ones(N,1)) , nonzero_signal_points(:,ones(N,1)).' ))*(floor(raster_width/2)*2+1);  
max_peak = cellfun(@(x) max(max(abs(x),[],2)), V).';  
minimum_peak = max_peak*max_peak.'.*sqrt_of_spikes*peak_discriminator_coeff; 

time_FFT = tic;
[aux_time1, corr1] = fast_corr(V{input_vec(1)},V{input_vec(2)}, T{input_vec(1)}, T{input_vec(2)}, dt, round(cor_delay_test/dt));
fprintf('time to evaluate FFTs = %.5g \n',toc(time_FFT))

figure(100), subplot(2,1,2),hold on,plot(-aux_time1, corr1)
plot(-aux_time1, minimum_peak(input_vec(1),input_vec(2))+0*aux_time1)

for j=1:length(sigma_test)
    iF_filtered = gaussfilt(aux_time1,corr1,sigma_test(j));
    hold on, plot(-aux_time1, iF_filtered),
end




