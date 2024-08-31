function estimated_delay = tdoa_excitationSource_method(mic1_aud, mic2_aud, index, windowsize, fs, shift)
    mic1_aud_lpr = LPres(mic1_aud,fs,20,5,10,1);
    mic2_aud_lpr = LPres(mic2_aud,fs,20,5,10,1);

    mic1_voiced_lpr = mic1_aud_lpr(index:index+windowsize-1);
    mic2_voiced_lpr = mic2_aud_lpr(index+shift: index+shift+windowsize-1);
    mic1_voiced_lpr_henv = HilbertEnv(mic1_voiced_lpr,fs);
    mic2_voiced_lpr_henv = HilbertEnv(mic2_voiced_lpr, fs);

    [estimated_delay, cc, peak_position, peak_value] = delay_gcc(mic1_voiced_lpr_henv, mic2_voiced_lpr_henv, "phat");


end