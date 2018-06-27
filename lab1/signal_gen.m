%function y_sig = signal_gen(sig_err_rate)
%y_sig = 0;
% Params:
USE_WARPLAB_TXRX        = 0;           % Enable WARPLab-in-the-loop (otherwise sim-only)
WRITE_PNG_FILES         = 0;           % Enable writing plots to PNG
CHANNEL                 = 11;          % Channel to tune Tx and Rx radios

% Waveform params
N_OFDM_SYMS             = 500;         % Number of OFDM symbols
MOD_ORDER               = 2;           % Modulation order (2/4/16/64 = BSPK/QPSK/16-QAM/64-QAM)
TX_SCALE                = 1.0;          % Scale for Tx waveform ([0:1])

% OFDM params
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16;                                     % Cyclic prefix length
N_DATA_SYMS             = N_OFDM_SYMS * length(SC_IND_DATA);      % Number of data symbols (one per data-bearing subcarrier per OFDM symbol)
INTERP_RATE             = 2;                                      % Interpolation rate (must be 2)

% WARPLab experiment params
USE_AGC                 = true;        % Use the AGC if running on WARP hardware
MAX_TX_LEN              = 2^20;        % Maximum number of samples to use for this experiment
TRIGGER_OFFSET_TOL_NS   = 3000;        % Trigger time offset toleration between Tx and Rx that can be accomodated

maximum_buffer_len = min(MAX_TX_LEN, 2^20);
    SAMP_FREQ           = 40e6;
    example_mode_string = 'sim';
%% Define a half-band 2x interpolation filter response remove
interp_filt2 = zeros(1,43);
interp_filt2([1 3 5 7 9 11 13 15 17 19 21]) = [12 -32 72 -140 252 -422 682 -1086 1778 -3284 10364];
interp_filt2([23 25 27 29 31 33 35 37 39 41 43]) = interp_filt2(fliplr([1 3 5 7 9 11 13 15 17 19 21]));
interp_filt2(22) = 16384;
interp_filt2 = interp_filt2./max(abs(interp_filt2));

% Define the preamble
% Note: The STS symbols in the preamble meet the requirements needed by the
% AGC core at the receiver. Details on the operation of the AGC are
% available on the wiki: http://warpproject.org/trac/wiki/WARPLab/AGC
sts_f = zeros(1,64);
sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
sts_t = ifft(sqrt(13/6).*sts_f, 64);
sts_t = sts_t(1:16);

% LTS for CFO and channel estimation
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);

% Use 30 copies of the 16-sample STS for extra AGC settling margin
preamble = [repmat(sts_t, 1, 30)  lts_t(33:64) lts_t lts_t];

% Sanity check variables that affect the number of Tx samples
num_samps_needed = ceil((TRIGGER_OFFSET_TOL_NS*1e-9) / (1/SAMP_FREQ)) + ...
                    INTERP_RATE*((N_OFDM_SYMS * (N_SC + CP_LEN)) + length(preamble) +  ceil(length(interp_filt2)/2));

                                
if(num_samps_needed > maximum_buffer_len)
    fprintf('Too many OFDM symbols for TX_NUM_SAMPS!\n');
    fprintf('Raise MAX_TX_LEN to %d, or \n', num_samps_needed);
    fprintf('Reduce N_OFDM_SYMS to %d\n', floor(((maximum_buffer_len - ceil((TRIGGER_OFFSET_TOL_NS*1e-9) / (1/SAMP_FREQ)))/INTERP_RATE - (length(preamble) +  ceil(length(interp_filt2)/2)))/(N_SC + CP_LEN)));
    return;
end
%% Generate a payload of random integers

%store the digital bits to tx_data.bin
tx_data = randi(MOD_ORDER, 1, N_DATA_SYMS) - 1;
complex_array = zeros(1, length(tx_data));
for j=1:length(tx_data)
    complex_array(j) = tx_data(j);
end    

fid = fopen('tx_data.bin', 'w');
fwrite(fid, complex_array, 'float');
fclose(fid);

% Functions for data -> complex symbol mapping (like qammod, avoids comm toolbox requirement)
% These anonymous functions implement the modulation mapping from IEEE 802.11-2012 Section 18.3.5.8
modvec_bpsk   =  (1/1)  .* [-1 1];
modvec_16qam  =  (1/sqrt(10)) .* [-3 -1 +3 +1];
modvec_64qam  =  (1/sqrt(43)) .* [-7 -5 -1 -3 +7 +5 +1 +3];

mod_fcn_bpsk  = @(x) complex(modvec_bpsk(1+x),0);
mod_fcn_qpsk  = @(x) complex(modvec_bpsk(1+bitshift(x, -1)), modvec_bpsk(1+mod(x, 2)));
mod_fcn_16qam = @(x) complex(modvec_16qam(1+bitshift(x, -2)), modvec_16qam(1+mod(x,4)));
mod_fcn_64qam = @(x) complex(modvec_64qam(1+bitshift(x, -3)), modvec_64qam(1+mod(x,8)));

% Map the data values on to complex symbols
switch MOD_ORDER
    case 2         % BPSK
        tx_syms = arrayfun(mod_fcn_bpsk, tx_data);
    case 4         % QPSK
        tx_syms = arrayfun(mod_fcn_qpsk, tx_data);
    case 16        % 16-QAM
        tx_syms = arrayfun(mod_fcn_16qam, tx_data);
    case 64        % 64-QAM
        tx_syms = arrayfun(mod_fcn_64qam, tx_data);
    otherwise
        fprintf('Invalid MOD_ORDER (%d)!  Must be in [2, 4, 16, 64]\n', MOD_ORDER);
        return;
end

% Reshape the symbol vector to a matrix with one column per OFDM symbol
tx_syms_mat = reshape(tx_syms, length(SC_IND_DATA), N_OFDM_SYMS);
%store the file tx_syms_mat.bin
complex_array_sym = zeros(1,length(tx_syms)*2);

for j=1:length(tx_syms)
    complex_array_sym(j) = real(tx_syms(j));
    complex_array_sym(j*2-1) = image(tx_syms(j));
end    

fid = fopen('tx_syms_mat.bin', 'w');
fwrite(fid, complex_array_sym, 'float');
fclose(fid);

% Define the pilot tone values as BPSK symbols
pilots = [1 1 -1 1].';

% Repeat the pilots across all OFDM symbols
pilots_mat = repmat(pilots, 1, N_OFDM_SYMS);
%% IFFT

% Construct the IFFT input matrix
ifft_in_mat = zeros(N_SC, N_OFDM_SYMS);

% Insert the data and pilot values; other subcarriers will remain at 0
ifft_in_mat(SC_IND_DATA, :)   = tx_syms_mat;
ifft_in_mat(SC_IND_PILOTS, :) = pilots_mat;

%Perform the IFFT
tx_payload_mat = ifft(ifft_in_mat, N_SC, 1);

% Insert the cyclic prefix
if(CP_LEN > 0)
    tx_cp = tx_payload_mat((end-CP_LEN+1 : end), :);
    tx_payload_mat = [tx_cp; tx_payload_mat];
end

% Reshape to a vector
tx_payload_vec = reshape(tx_payload_mat, 1, numel(tx_payload_mat));

% Construct the full time-domain OFDM waveform
tx_vec = [preamble tx_payload_vec];

% Pad with zeros for transmission to deal with delay through the
% interpolation filter
tx_vec_padded = [tx_vec, zeros(1, ceil(length(interp_filt2)/2))];
%% Interpolate
% Zero pad then filter (same as interp or upfirdn without signal processing toolbox)
if( INTERP_RATE ~= 2)
   fprintf('Error: INTERP_RATE must equal 2\n'); 
   return;
end

tx_vec_2x = zeros(1, 2*numel(tx_vec_padded));
tx_vec_2x(1:2:end) = tx_vec_padded;
tx_vec_air = filter(interp_filt2, 1, tx_vec_2x);


% Scale the Tx vector to +/- 1: Normalization
tx_vec_air = TX_SCALE .* tx_vec_air ./ max(abs(tx_vec_air));

TX_NUM_SAMPS = length(tx_vec_air);

if(USE_WARPLAB_TXRX)
    wl_basebandCmd(nodes, 'tx_delay', 0);
    wl_basebandCmd(nodes, 'tx_length', TX_NUM_SAMPS);                                                        % Number of samples to send
    wl_basebandCmd(nodes, 'rx_length', TX_NUM_SAMPS + ceil((TRIGGER_OFFSET_TOL_NS*1e-9) / (1/SAMP_FREQ)));   % Number of samples to receive
end
%% Signal with Noise

rx_vec_air = [tx_vec_air, zeros(1,ceil((TRIGGER_OFFSET_TOL_NS*1e-9) / (1/SAMP_FREQ)))];
 % 5, 10, 15, 20, 25
SNR = 5;
Pn = mean(abs(rx_vec_air).^2)/SNR;
rx_vec_air = rx_vec_air + sqrt(Pn) * complex(randn(1,length(rx_vec_air)), randn(1,length(rx_vec_air))) / sqrt(2);

% store the rx_vec_air.bin
complex_array_rx = zeros(1, 2*length(rx_vec_air));
for j=1:length(rx_vec_air)
    complex_array_rx(2*j-1) = real(rx_vec_air(j));
    complex_array_rx(2*j) = imag(rx_vec_air(j));
end    

fid = fopen('rx_vec_air.bin', 'w');
fwrite(fid, complex_array_rx, 'float');
fclose(fid);

%end