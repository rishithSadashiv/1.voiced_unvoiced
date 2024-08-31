function tdoa = samples_to_time(samples, fs)
    tau = samples / fs;
    tdoa = tau*1000;

end