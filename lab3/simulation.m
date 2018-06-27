%% call the function signal_gen and decode for each SNR
clear;
close all;
clc;

signal_gen(5);
[ snr, snr1, snr2, snr3, sig1 ] = decode();
y1 = [snr snr1 snr2 snr3];

signal_gen(10);
[ snr, snr1, snr2, snr3, sig2 ] = decode();
y2 = [snr snr1 snr2 snr3];

signal_gen(15);
[ snr, snr1, snr2, snr3, sig3 ] = decode();
y3 = [snr snr1 snr2 snr3];

signal_gen(20);
[ snr, snr1, snr2, snr3, sig4 ] = decode();
y4 = [snr snr1 snr2 snr3];

signal_gen(25);
[ snr, snr1, snr2, snr3, sig5 ] = decode();
y5 = [snr snr1 snr2 snr3];

x = 5:5:25;
y = [y1; y2; y3; y4; y5];

figure;
plot(abs(sig1));
title("Signal Amplitude (SNR 5)");

figure;
plot(abs(sig2));
title("Signal Amplitude (SNR 10)");

figure;
plot(abs(sig3));
title("Signal Amplitude (SNR 15)");

figure;
plot(abs(sig4));
title("Signal Amplitude (SNR 20)");

figure;
plot(abs(sig5));
title("Signal Amplitude (SNR 25)");

figure;
bar(x,y);
grid on;
title('SNR Figure');
xlabel('Simulated SNR (dB)');
ylabel('Actual SNR');

