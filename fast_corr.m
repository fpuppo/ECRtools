%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020, Francesca Puppo
%
% This function computes cross-correlations 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [aux_time, corr] = fast_corr(V1,V2, T1, T2, dt, size_delay)

data1tt = round(T1 / dt); 
data2tt = round(T2 / dt); 

T2 = data2tt(1,:);  

aux_ind1 = [ones(size(data1tt)), 2 * ones(size(data2tt))];
aux_ind2 = [data1tt, data2tt];
aux_val = [V1, V2]; 
spar_rec = sparse(aux_ind1, aux_ind2, aux_val,2,max(max(aux_ind2))+size_delay); 

num_ind = size(T2, 2); 

aux_time = (- size_delay: size_delay); 
aux_mat = ones(num_ind, 1) * aux_time; 
aux_mat = aux_mat + T2(1, :).'; 

test_spar = spar_rec(1, :); 
aux_signal = test_spar(aux_mat); 

vv_extend = V2; vv_extend(size(aux_signal, 2), 1) = 0; 

fft_signal = fft(full(aux_signal)'); 
fft_vv = fft(vv_extend); 
fft_prod = conj(fft_signal) .* fft_vv; 
corr = sum(ifft(fft_prod),2); 

aux_time = -aux_time*dt+dt;
aux_time = aux_time(end:-1:1).';
corr = corr(end:-1:1);


