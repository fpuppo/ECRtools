classdef Spike_v1 < Waveform
    %SPIKE_V1: An extension of Waveform that represents a spike recorded by
    %          spike detector in Axis.
    %
    %   TriggerSampleOffset: Offset (in samples) from the start of the
    %                        waveform where the spike detector was
    %                        triggered
    %
    %   StandardDeviation:   RMS voltage value of the signal noise at the
    %                        time the spike was caputred
    %
    %   ThresholdMultiplier: Multiplier(if applicable) of the RMS Noise that was 
    %                        used as the trigger voltage for this spike
    %
    
    
    properties(GetAccess = public, Constant = true)
        LOADED_HEADER_SIZE = 30;
    end
    
    properties (GetAccess = public, SetAccess = private)
        TriggerSampleOffset
        StandardDeviation
        ThresholdMultiplier
    end
    
    methods
        
        function this = Spike_v1( ...
                aChannel, ...
                aStart, ...
                aData, ...
                aSource,  ...
                aTriggerSampleOffset, ...
                aStandardDeviation, ...
                aThresholdMultiplier)
            if(nargin == 0)
                return;
            end
            
            this.Channel = aChannel;
            this.Start = aStart;
            this.Data = aData;
            this.Source = aSource;
            this.TriggerSampleOffset = aTriggerSampleOffset;
            this.StandardDeviation = aStandardDeviation;
            this.ThresholdMultiplier = aThresholdMultiplier;
            %We do not support spike data subsampling at this time.
            this.SubsampleFactor = 1;
        end
    end
    
end

