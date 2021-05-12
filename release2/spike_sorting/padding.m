function [spk_wave, spk_wave_raster, time_pad] = padding(neuron_spk_V, neuron_spk_T, maxt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LAST UPDATE 06/29/2020, by Francesca Puppo
%
%%%%% Pre-processing: 
% individual spikes are concatenated (zeros in between spikes) to reconstruct temporal sequence of spike events;
% neurons' activities are aligned in time and lhs and rhs are padded with
% zeros in preparation to Fourier analysis (avoid border effects)
%
% Outputs:
% time_pad: reference series of time stamps
% spk_wave: single spike waveforms ad zeros in between
% spk_wave_raster: unit steps correspond to single spikes; zeros in between
% spikes

spk_wave = cell(1,size(neuron_spk_T,1));
spk_wave_raster = cell(1,size(neuron_spk_T,1));
t_wave = cell(1,size(neuron_spk_T,1));


for neu=1:size(neuron_spk_T,1)
    if size(neuron_spk_T{neu},2)
    % for each neuron, concatenate all spikes, padding with zeros in between spikes 
    for spknum = 1:size(neuron_spk_T{neu},2)
               
        dt = neuron_spk_T{neu}(2,spknum) - neuron_spk_T{neu}(1,spknum);
        
        s = [0; neuron_spk_V{neu}(:,spknum); 0];
        s_raster = abs(s)== max(abs(s));
        t = [neuron_spk_T{neu}(1,spknum)-dt; neuron_spk_T{neu}(:,spknum); neuron_spk_T{neu}(end,spknum)+dt];        
        
        spk_wave{neu} = [spk_wave{neu}; s];
        spk_wave_raster{neu} = [spk_wave_raster{neu}; s_raster];
        t_wave{neu} = [t_wave{neu}; t];
        
    end
    
    % construct a reference time vector based on the maximum spike time
    % detected across all elctrodes in a well (maxt)
    time_pad = linspace(-maxt/2,maxt+(maxt/2),2*round((maxt)/dt)+1).';     % pad with zeros the lhs and rhs to avoid border effects in FFT

    [t_wave{neu},I] = sort(t_wave{neu});
    spk_wave{neu} = spk_wave{neu}(I);
    spk_wave_raster{neu} = spk_wave_raster{neu}(I);
    [t_wave{neu}, I] = unique(t_wave{neu},'last');
    spk_wave{neu} = spk_wave{neu}(I);
    spk_wave_raster{neu} = spk_wave_raster{neu}(I);
    
    spk_wave{neu} = interp1([time_pad(1); t_wave{neu}; time_pad(end)],[0; spk_wave{neu}; 0 ],time_pad);
%     test = sum(spk_wave_raster{neu});
    spk_wave_raster{neu} = interp1([time_pad(1); t_wave{neu}; time_pad(end)],[0; spk_wave_raster{neu}; 0 ],time_pad);
    spk_wave_raster{neu} = spk_wave_raster{neu}>.4;
%     if sum(spk_wave_raster{neu})~=test
%         [sum(spk_wave_raster{neu}) test]
%         error('raster failure')
%     end
    end
end

