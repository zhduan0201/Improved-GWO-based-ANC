function [y,noise] = add_awgn(x,snr) %x: original signal, y: noisy signal, snr: signal-to-noise ratio
noise=randn(size(x));
Nx=length(x);
signal_power=1/Nx*sum(x.*x);
noise_power=1/Nx*sum(noise.*noise);
noise_variance=signal_power/(10^(snr/10));
noise=sqrt(noise_variance/noise_power)*noise;
y=x+noise;
