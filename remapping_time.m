%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020, Francesca Puppo
%
% This function loads spike sorted data, checks data consistency, aligns
% recordings in time and pad with zeros if not previously done
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ref_time, V, T] = remapping_time(folder,h,raster_width,cor_delay_interval)

%% -----------------------------------------------------------------
% read contents in folder and load spike sorted data contained in files that have file name including the characters "WELL"
file = dir(folder);
[name{1:length(file)}] = deal(file.name);
name = name(cellfun(@(x) contains(x,['raster_spkwave_WELL' num2str(h) '_']),name));

%% -----------------------------------------------------------------
% load data and check consistency
max_time = -inf;
min_time = inf;
dt = zeros(length(name),1);
for j = 1:length(name)
    load([folder name{j}],'time_pad')
    if length(time_pad) ==3 ; time_pad = linspace(time_pad(1),time_pad(2),time_pad(3)); end
    if ~isempty(time_pad)
        if time_pad(end)>max_time
            max_time = time_pad(end);
        end
        if time_pad(1)<min_time
            min_time = time_pad(1);
        end
        dt(j) = mean(diff(time_pad));
        % check data consistency: uniform time discretization along the entire recording
        if any(abs(dt(j) - diff(time_pad))>1e-10/(max_time-min_time)*length(time_pad))
            error('time discretization not uniform')
        end
    else
        dt(j) = nan;
    end
end

%% -----------------------------------------------------------------
% check data consistency: the activities of different neurons must have same time discretization
if all(isnan(dt)), error('no activity data loaded'), end
if any( abs( mean(dt(~isnan(dt)))-dt(~isnan(dt)) ) > 1e-10/(max_time-min_time)*length(time_pad))
    error('time discretization not cosistent among different neurons')
end
dt = dt(~isnan(dt));
dt = dt(1);

%% -----------------------------------------------------------------
% Align activities in time and prepare data for Fourier analysis
% (all recordings must be padded with zeros on each side to avoid border effects with calculation of fft; all must have the same length)

% define the reference time
ref_time = (min_time:dt:max_time+dt/2).'; % dt/2 is included in the rhs to make sure max_time is included as last point

% pad  neurons' spiking activities with zeros
k = 0;
V = {};
T = {};
delta_t = ceil(max(cor_delay_interval)*2/dt);
for j = 1:length(name)
    load([folder name{j}],'time_pad','spk_wave_raster')
    spk_wave = spk_wave_raster;
    if length(time_pad) ==3 ; time_pad = linspace(time_pad(1),time_pad(2),time_pad(3)); end
    time_pad = [time_pad(:); time_pad(end)+(time_pad(2)-time_pad(1))*(1:delta_t).'].';
    for h=1:length(spk_wave)
        ind = find(spk_wave{h}>.5);
        ind_vec = delta_t+ceil(raster_width/2)+[ind+(-floor(raster_width/2):-1) ind ind+(1:floor(raster_width/2))];
        ind = find(sparse(1,ind_vec,1)==2);
        V{end+1} = ones(size(ind_vec)).';
        for jj=1:size(ind_vec,2)
            V{end}( jj, builtin('_ismemberhelper',ind_vec(:,jj),ind)) = 1/2;
        end
        T{end+1} = time_pad(ind_vec).';
    end
    
    
end




