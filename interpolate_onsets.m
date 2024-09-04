function [sine_out, saw_out] = interpolate_onsets (onsets) 

% This function imports onsets from the "drifting metronomes experiment"
% and returns the interpolated timeseries in the form of a sinewave (for
% e.g. JRQA) and rampwave (i.e., 'chain-saw' phase timeseries)
%
% NOTES:
% - The sampling rate for interpolation must be 1000Hz, because onsets have
%   1ms time resolution. Any other sampling rate results in approximation
%   and misalagnment. Downsampling can be performed on the output, AFTER
%   this interpolation process.
% - This approach sacrifices the very first and the very last tapping onsets;
%   the reason is that subjects will never tap exactly at the same time, so
%   we cannot take their behaviour as reference for defining the timeseries
%   boundaries. This would inevitably result in misalignment.
%
% Guaranteing the timeseries alignment is absolute priority for reliable
% results.

        % Recording settings        
        srate = 1000; % time resolution from Teensy = 1ms
        time  = 390;    % total duration, in seconds
        
        % Pre-allocate output (create a 'grid' of 0s, to fill with taps
        % positions; size in known beforehand and equal for every participant)
        saw_out = zeros(time*srate, 1); 
        
        % Put taps in position on the 0s grid, for interpolation
        %saw_out(1:onsets(1)) = nan;
        saw_out(onsets(2:end),1) = 1; %skip first tap, which might be negative (hence outside the vector)
        idx = find(saw_out~=0);       % get index of onsets in the grid for interpolation 
        
        % NaN padding (we use non-numerical values, not to bias the results) 
        saw_out(1:idx(1)-1) = nan;         % sacrifice samples up to first tap
        saw_out(idx(end):end) = nan;    % sacrifice samples from last tap until the end
       
       
        % Compute time series
        for i = 1:length(idx)-1
            
            % Linearly interpolate segments with ramp waves from 0 to 2pi
            saw_out(idx(i):idx(i+1)) = linspace( 0 , 2*pi , length(idx(i):idx(i+1)));
            
        end
        
        % Compute sinewave
        sine_out = sin(saw_out);
        

end