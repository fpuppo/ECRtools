function [neuron_spk_V, neuron_spk_T] = cluster_data(t,v,min_detected_spikes,exvariance_mode,cluster_plot,sorted_spike_plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LAST UPDATE 06/11/2020, by Francesca Puppo
%
%%%%% Spike sorting based on PCA and k-means clustering 
%

%% PCA analysis

% discriminant = zeros(size(v,2),1);
% v1 = v;
% tt = (1:size(v,1)).';
% for j=1:size(v1,2)
%     coeff = [tt ones(size(v1,1),1)]\v1(:,j);
%     v1(:,j) = v1(:,j)-coeff(1)*tt-coeff(2);
%     discriminant(j) = std(v1(:,j)-coeff(1)*tt-coeff(2));
% end
[~, score, latent] = pca(v.');%,'Centered',false,'Weights',discriminant(discriminant>max(discriminant)/5));
% score = v1.'/(COEFF.');

D =delaunayTriangulation(score(:,1:3));
p1 = D.Points(:,1);
p2 = D.Points(:,2);
p3 = D.Points(:,3);
pairs = unique(sort([D.ConnectivityList(:,[1 2]); D.ConnectivityList(:,[2 3]); D.ConnectivityList(:,[3 4]); D.ConnectivityList(:,[4 1])],2),'rows');
distances = sqrt(diff(p1(pairs),[],2).^2+diff(p2(pairs),[],2).^2+diff(p3(pairs),[],2).^2);
points = pairs(:);
distances = [distances; distances];

% minimum_distance = full(1./max(sparse(1:length(points),points,1./distances),[],1));
average_distance = full(sum(sparse(1:length(points),points,distances),1)./sum(sparse(1:length(points),points,1),1));
th = 2*mean(average_distance)-min(average_distance);
v = v(:,average_distance<th) ;
%v1 = v1(:,average_distance<th) ;
t = t(:,average_distance<th) ;
score = score(average_distance<th,:) ;

if exvariance_mode == 1
    figure(3)
    plot(latent)
    title('explained variance')
end

%% k-means clustering

% NOTE: Data reported in the preprint article came from
% MEA with fewer electrodes and therefore less spikes to
% analyze. Clustering was performed with the k-means
% approach and k was visually inspected.
% Cleber's data were recorded on higher density MEA, therefore resulting
% in many more neurons and more spikes to analyze.
% For this dataset, evaluating k manually is practically impossible.
% Therefore I decided to use the gap criterion to find
% k.
% The gap criterion has been used for spike sorting in
% previous works by other groups.
% However, different clustering approaches could be
% investigated and might perform better. More in general, different
% spike-sorting algorithms could be also investigated to improve
% computation efficiency and clustering accuracy

% Gap evaluation criterion for k

lastwarn('')
n_score = min(7,size(v,2)-1);    
eva = evalclusters(score(:, 1:n_score),'kmeans','gap','KList',1:min(n_score,6),'ReferenceDistribution','PCA');
k = eva.OptimalK;
[~, LASTID] = lastwarn;
if ~isempty(lastwarn)&&(~strcmp(LASTID,'stats:kmeans:FailedToConvergeRep'))
    error('vedi che il gap criterion e'' un cretino')
end
% clustering
if k>1
    lastwarn('')
    [cluster_idx, ~, ~] = kmeans(score(:, 1:n_score), k,'Replicates', 100,'MaxIter',1000);%,'Distance','correlation' );
    if ~isempty(lastwarn)
        warning('vedi che kmeans e'' un cretino')
    end
    
    %[s,h] = silhouette(score(:, 1:n_score),cluster_idx);
    
    color = ['k' 'g' 'r' 'b' 'm' 'c' 'y'];
    neuron_spk_T = cell(k,1);
    neuron_spk_V = cell(k,1);
    
    spk_ind = cell(k,1);
    for i = 1:k
        matching_points = cluster_idx == i;
        if length(matching_points)>min_detected_spikes         % discard neurons with less than a minimum number of spikes that the user can select (if less, the neuron is considered inactive)
            spk_ind{i} = find(cluster_idx==i);
            neuron_spk_T{i} = t(:,spk_ind{i});
            neuron_spk_V{i} = v(:,spk_ind{i});
        else
            neuron_spk_T{i} = [];
            neuron_spk_V{i} = [];
        end
    end
else
    matching_points = 1:size(v,2);
    cluster_idx = 1;
    neuron_spk_T{k} = t;
    neuron_spk_V{k} = v;
end

if cluster_plot==1
    figure(4)
    for i =1:k
        matching_points = cluster_idx == i;
        scatter3(score(matching_points, 1), score(matching_points, 2), score(matching_points, 4),color(i))
        hold on
    end
end

% bad implementation; just to visualize spike sorted data
if sorted_spike_plot==1
    figure(6)
    for j=1:k
        subplot(k,1,j), plot(neuron_spk_V{j},color(j))
    end    
    xlabel('time (s)');ylabel('voltage (V)');
%     title(['WELL' num2str(well_R(wr)) ',' num2str(well_C(wc)) '- Electrode' num2str(eltr_C(elc)) ',' num2str(eltr_R(elr))])
end

