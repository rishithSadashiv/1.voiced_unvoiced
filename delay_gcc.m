function [estimated_delay,cc,peak_position,peak_value] = delay_gcc( input1, input2, method )
    % Code taken from https://github.com/javiribera/TDE-and-whale-localization/tree/master/algorithms_TDE
    % DELAY_GCC Estimates the delay of 2 signals using the function GCC.m,
    % the weighting filter METHOD and taking the maximum
    
    if length(input1) ~= length(input2)
        error('- ¡Both inputs must be the same length to correlate them!');
    end
    cc = gcc(input1,input2, method);
    cc = cc/max(cc);
    [value, peak_positions] = max(cc);
    peak_position = max(peak_positions);
    peak_value = value;
    estimated_delay = peak_position - length(input1);
end
