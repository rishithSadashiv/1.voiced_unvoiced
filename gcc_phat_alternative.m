% gcc_phat other implementation
function [estimated_delay,GCCn, peak_position] = gcc_phat_alternative(mic1_segment, mic2_segment)
    wlen = length(mic1_segment);
    w = hann(wlen);
    mic1_voiced = mic1_segment.*w;
    mic2_voiced = mic2_segment.*w;
    
    Nfft = wlen;
    GCC = fftshift(real(ifft(exp(1i*angle(fft(mic1_voiced,Nfft).*conj(fft(mic2_voiced,Nfft))))))); 
    
    GCCn = GCC/max(GCC);
    
    [value, peak_positions] = max(GCC);
    peak_position = max(peak_positions);
    estimated_delay = (peak_position - wlen/2)-1;


end
