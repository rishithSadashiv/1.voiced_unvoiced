% Both audio files are recorded using ambisonic microphone in Ambisonic A
% mode. The distance was 15cm. And azimuth angle is 30degrees


clc;
clear all;
close all;

[wav1,fs]=audioread('./20240901/15cm/30degree/mic1_15cm_30deg.WAV');
[wav2,fs]=audioread('./20240901/15cm/30degree/mic2_15cm_30deg.WAV');

wav1=wav1(:,2);
wav2=wav2(:,4);

wav1=wav1-mean(wav1);
wav1=wav1./max(abs(wav1));
wav2 = wav2-mean(wav2);
wav2=wav2./max(abs(wav2));

ms30 = 0.03 * fs;
voiced_index = 25.00*fs;
unvoiced_index = 24.320*fs;

%%%%%%%%%%%%%%%% Excitation source method

mic1_aud_lpr = LPres(wav1,fs,20,5,10,1);
mic2_aud_lpr = LPres(wav2,fs,20,5,10,1);


mic1_voiced_lpr = mic1_aud_lpr(voiced_index:voiced_index+ms30-1);
mic2_voiced_lpr = mic2_aud_lpr(voiced_index: voiced_index+ms30-1);
mic1_voiced_lpr_henv = HilbertEnv(mic1_voiced_lpr,fs);
mic2_voiced_lpr_henv = HilbertEnv(mic2_voiced_lpr, fs);

disp("Using Excitation Information");
[estimated_delay, cc, peak_position, peak_value] = delay_gcc(mic1_voiced_lpr_henv, mic2_voiced_lpr_henv, "cc");
disp("Voiced segment: estimated delay = " + estimated_delay)

mic1_unvoiced_lpr = mic1_aud_lpr(unvoiced_index:unvoiced_index+ms30-1);
mic2_unvoiced_lpr = mic2_aud_lpr(unvoiced_index: unvoiced_index+ms30-1);
mic1_unvoiced_lpr_henv = HilbertEnv(mic1_unvoiced_lpr, fs);
mic2_unvoiced_lpr_henv = HilbertEnv(mic2_unvoiced_lpr, fs);

[estimated_delay, cc, peak_position, peak_value] = delay_gcc(mic1_unvoiced_lpr_henv, mic2_unvoiced_lpr_henv, "cc");
disp("Unvoiced segment: estimated delay = " + estimated_delay)


%%%%%%%%%%%%%%%%% GCC-PHAT

mic1_voiced = wav1(voiced_index:voiced_index+ms30-1);
mic2_voiced = wav2(voiced_index:voiced_index+ms30-1);

hamming_w = hamming(length(mic1_voiced));
mic1_voiced_w = mic1_voiced.*hamming_w;
mic2_voiced_w = mic2_voiced.*hamming_w;

mic1_unvoiced = wav1(unvoiced_index:unvoiced_index+ms30-1);
mic2_unvoiced = wav2(unvoiced_index:unvoiced_index+ms30-1);


disp("Using GCC-PHAT");
[estimated_delay, cc, peak_position, peak_value] = delay_gcc(mic1_voiced_w, mic2_voiced_w, "phat");
disp("Voiced segment: estimated delay = " + estimated_delay)

[estimated_delay, cc, peak_position, peak_value] = delay_gcc(mic1_unvoiced, mic2_unvoiced, "phat");
disp("Unvoiced segment: estimated delay = " + estimated_delay)












