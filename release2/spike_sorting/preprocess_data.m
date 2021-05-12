function [out_waveform_mat, out_time_mat] = preprocess_data(v,t)

% figure
% plot(v)

%% filter out 'flat' spike waveforms
negid = find(all(v<0,1));       % find waveforms with only negative values
posid = find(all(v>0,1));       % find waveforms with only positive values

if ~isempty(negid) || ~isempty(posid)
v(:,[negid posid]) = [];
t(:,[negid posid]) = [];
end

out_waveform_mat = v;
out_time_mat = t;

% 
% %% filter out low-prominence spike waveforms
% 
% max(max(v))
% min(min(v))
% 
% 
% [I,J] = find(v==min(v,[],1));
% [J,IA,~] = unique(J);
% I = I(IA);
% v1 = v(1:19,:);
% for j=1:size(v,2)
%     ind = I(j)+(-9:9);
%     m0 = ind>0;
%     m39 = ind<39;
%     ind = ind((m0)&(m39));
%     v1(1:19,j) = [zeros(sum(~m0),1); v(ind,j); zeros(sum(~m39),1)];
% end
% discriminant = zeros(size(v1,2),1);
% v2 = v1;
% tt = (1:19).';
% for j=1:size(v1,2)
%     coeff = [tt ones(size(v1,1),1)]\v1(:,j);
%     v2(:,j) = v1(:,j)-coeff(1)*tt-coeff(2);
%     discriminant(j) = std(v1(:,j)-coeff(1)*tt-coeff(2));
%     
% end
% [~, I] = sort(discriminant);
% [~, score, latent] = pca(v1(:,I(end:-1:round(0.3*length(I)))).');
% eva = evalclusters(score,'kmeans','gap','KList',[1:6],'ReferenceDistribution','PCA');
% 
% 
% n_score = 10;
% if n_score>=size(v1,1)
%     n_score = size(v1,1)-1;
% end
% [cluster_idx, ~, ~] = kmeans(score(:, 1:n_score), k,'Replicates', 100);
% 
% color = ['k' 'g' 'r' 'b' 'm' 'c' 'y'];
% neuron_spk_T = cell(k,1);
% neuron_spk_V = cell(k,1);
% 
% spk_ind = cell(k,1);
% for i = 1:k
%     matching_points = cluster_idx == i;
%     if length(matching_points)>min_detected_spikes         % discard neurons with less than a minimum number of spikes that the user can select (if less, the neuron is considered inactive)
%         spk_ind{i} = find(cluster_idx==i);
% %         neuron_spk_T{i} = t(:,spk_ind{i});
%         neuron_spk_V{i} = v1(spk_ind{i},:);
%     else
% %         neuron_spk_T{i} = [];
%         neuron_spk_V{i} = [];
%     end
% end
% 
% if cluster_plot==1
%     figure(4)
%     for i =1:k
%         matching_points = cluster_idx == i;
%         scatter3(score(matching_points, 1), score(matching_points, 2), score(matching_points, 3),color(i))
%         hold on
%     end
% end
% 
% % bad implementation; just to visualize spike sorted data
% 
%     figure(5)
%     for j=1:k
%         subplot(k,1,j), plot(neuron_spk_V{j},color(j))
%     end    
% 
% 
% 
% 
% 
% 
% % linear regression
% discriminant = zeros(size(v,2),1);
% v1 = v;
% for j=1:size(t,2)
%     coeff = [t(:,j) ones(size(t,1),1)]\v(:,j);
%     v1(:,j) = v(:,j)-coeff(1)*t(:,j)-coeff(2);
%     discriminant(j) = std(v(:,j)-coeff(1)*t(:,j)-coeff(2));
%     
% end
% [~, I] = sort(discriminant);
% [~, score, latent] = pca(v(:,I(end:-1:round(0.5*length(I)))).');
% eva = evalclusters(score,'kmeans','gap','KList',[1:6],'ReferenceDistribution','PCA');
% 
% 
% 
% p = running_percentile(v(:,1),100,10);
% 
% figure, plot(v(:,1)), hold on, plot(p)


