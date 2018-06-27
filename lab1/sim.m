%% call the function signal_gen and decode for each SNR
%{
clear;
close all;
clc;

SNR = [5 10 15 20 25];
stat_table = [];

for j=1:5
   signal_gen(SNR(j));
   [ber,snr] = decode();
   stat_table = [stat_table; ber snr];
end

x = 5:5:25;
y_ber = stat_table(:,1);
y_snr = stat_table(:,2);

figure;
bar(x,y_ber,'r');
grid on;
title('BER Figure');
xlabel('Simulated SNR (dB)');
ylabel('BER');

figure;
bar(x,y_snr,0.5);
grid on;
title('SNR Figure');
xlabel('Simulated SNR (dB)');
ylabel('Actual SNR');
%}
rx_vec_air_mix = read_complex_binary( 'rx_signals.bin' );
figure;
plot(angle(rx_vec_air_mix));
