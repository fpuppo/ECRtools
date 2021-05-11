%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020, Francesca Puppo
%
% This function performs the 'super-selective' filtering of apparent and
% indirect connections
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [adjmat, delay_adjmat] = edge_covering_minimization(delays, peaks, epsilon)

adjmat = zeros(size(peaks));
delay_adjmat = zeros(size(peaks));
for j = 1:size(delays,1)                                                   % go through all neurons
    for k =setdiff(1:size(delays,2),j)                                     % exclude autocorrelations
        direct_correlation_index = false(size(delays{j,k}));
        for hjk = 1:length(delays{j,k})                                    % index of correlation peaks between neurons j and k
            direct_correlation_index(hjk) = check_triangle(delays,epsilon,peaks,j,k,hjk);   
        end
        filtered_delays = delays{j,k}(direct_correlation_index);
        filtered_delays = filtered_delays(filtered_delays<0);
        adjmat(j,k) = ~isempty(filtered_delays);
        delay_adjmat(j,k) = min([-filtered_delays; inf]);
    end
end


function direct_correlation_index = check_triangle(delays,epsilon,peaks,j,k,hjk)

direct_correlation_index = true;
for m = 1:size(delays,2)                                       % index of the third neuron to close the triangle
    if (m ~= j) && (m ~= k)                                    % third neuron must be different from j and k neurons
        for hjm = 1:length(delays{j,m})                        % index of correlation peaks between neurons j and m
            for hmk = 1:length(delays{m,k})                    % index of correlation peaks between neurons m and k
               if abs((delays{j,m}(hjm) + delays{m,k}(hmk) - delays{j,k}(hjk))) <= epsilon  % check triangle with cyclec indexing
                    Pvec = [peaks{j,m}(hjm) peaks{m,k}(hmk) peaks{j,k}(hjk)];
                    minid = find(min(Pvec)==Pvec);  %--> remove connection
                    if minid==3
                        direct_correlation_index = false;
                        return
                    end
                end
            end
        end
    end
end



