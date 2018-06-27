SFO = decode(1,1);

non_SFO = decode(0,0);

fi_1 = fopen('tx_vec.bin','rb');
x_inter_1 = fread(fi_1, 'float32');

fi_2 = fopen('rx_signals.bin','rb');

x_inter_2 = fread(fi_2, 'float32');

% if data is complex
x_1 = x_inter_1(1:2:end) + 1i*x_inter_1(2:2:end);
figure;
plot(abs(x_1).^2);
title("Tx Data");

rx_vec_air = read_complex_binary('rx_signals.bin');
offset_tx = 108400;
raw_rx_dec = rx_vec_air(offset_tx:12000+offset_tx).';


x_2 = x_inter_2(1:2:end) + 1i*x_inter_2(2:2:end);
figure;
plot(abs(x_2).^2);
title("Raw Rx Signal");

figure;
plot(abs(raw_rx_dec).^2);
title("Filtered Rx Signal");

figure; 

plot(non_SFO, 'x');
hold on;
plot(SFO, 'o');
ylim([0 3.5]);
title('Phases of decoded signal');


