function [ output ] = signal_gen( SNR_simulated )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % wl_example_siso_ofdm_txrx.m
    % A detailed write-up of this example is available on the wiki:
    % http://warpproject.org/trac/wiki/WARPLab/Examples/OFDM
    %
    % Copyright (c) 2015 Mango Communications - All Rights Reserved
    % Distributed under the WARP License (http://warpproject.org/license)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Tx side
    % Params:
    % Change into simualtion mode.
    USE_WARPLAB_TXRX        = 0;           % Enable WARPLab-in-the-loop (otherwise sim-only)
    CHANNEL                 = 11;          % Channel to tune Tx and Rx radios

    % Waveform params
    N_OFDM_SYMS             = 50;         % Number of OFDM symbols
    MOD_ORDER               = 2;           % Modulation order (2/4/16/64 = BSPK/QPSK/16-QAM/64-QAM)
    TX_SCALE                = 1.0;         % Scale for Tx waveform ([0:1]) (TX gain)

    % OFDM params
    SC_IND_PILOTS           = [ 8 , 22 , 44 , 58 ];                                             % Pilot subcarrier indices
    % Index 1 is for DC and index 28 to 28 is not usesd form leakage.
    SC_IND_DATA             = [ 2 : 7 , 9 : 21 , 23 : 27 , 39 : 43 , 45 : 57 , 59 : 64 ];       % Data subcarrier indices
    N_SC                    = 64;                                                               % Number of subcarriers
    CP_LEN                  = 16;                                                               % Cyclic prefix length
    N_DATA_SYMS             = N_OFDM_SYMS * length(SC_IND_DATA);                                % Number of data symbols (one per data-bearing subcarrier per OFDM symbol)
    INTERP_RATE             = 2;                                                                % Interpolation rate (must be 2)

    % WARPLab experiment params
    USE_AGC                 = true;        % Use the AGC if running on WARP hardware
    MAX_TX_LEN              = 2 ^ 20;      % Maximum number of samples to use for this experiment
    TRIGGER_OFFSET_TOL_NS   = 3000;        % Trigger time offset toleration between Tx and Rx that can be accomodated


    if USE_WARPLAB_TXRX

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set up the WARPLab experiment
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        NUMNODES = 2;

        % Create a vector of node objects
        nodes   = wl_initNodes(NUMNODES);
        node_tx = nodes(1);
        node_rx = nodes(2);

        % Create a UDP broadcast trigger and tell each node to be ready for it
        eth_trig = wl_trigger_eth_udp_broadcast;
        wl_triggerManagerCmd(nodes, 'add_ethernet_trigger', [eth_trig]);

        % Read Trigger IDs into workspace
        trig_in_ids  = wl_getTriggerInputIDs(nodes(1));
        trig_out_ids = wl_getTriggerOutputIDs(nodes(1));

        % For both nodes, we will allow Ethernet to trigger the buffer baseband and the AGC
        wl_triggerManagerCmd(nodes, 'output_config_input_selection', [trig_out_ids.BASEBAND, trig_out_ids.AGC], [trig_in_ids.ETH_A]);

        % Set the trigger output delays.
        nodes.wl_triggerManagerCmd('output_config_delay', [trig_out_ids.BASEBAND], 0);
        nodes.wl_triggerManagerCmd('output_config_delay', [trig_out_ids.AGC], TRIGGER_OFFSET_TOL_NS);

        % Get IDs for the interfaces on the boards. 
        ifc_ids_TX = wl_getInterfaceIDs(node_tx);
        ifc_ids_RX = wl_getInterfaceIDs(node_rx);

        % Set up the TX / RX nodes and RF interfaces
        TX_RF     = ifc_ids_TX.RF_A;
        TX_RF_VEC = ifc_ids_TX.RF_A;
        TX_RF_ALL = ifc_ids_TX.RF_ALL;

        RX_RF     = ifc_ids_RX.RF_A;
        RX_RF_VEC = ifc_ids_RX.RF_A;
        RX_RF_ALL = ifc_ids_RX.RF_ALL;

        % Set up the interface for the experiment
        wl_interfaceCmd(node_tx, TX_RF_ALL, 'channel', 2.4, CHANNEL);
        wl_interfaceCmd(node_rx, RX_RF_ALL, 'channel', 2.4, CHANNEL);

        wl_interfaceCmd(node_tx, TX_RF_ALL, 'tx_gains', 3, 30);

        if USE_AGC
            wl_interfaceCmd(node_rx, RX_RF_ALL, 'rx_gain_mode', 'automatic');
            wl_basebandCmd(nodes, 'agc_target', -13);
        else
            wl_interfaceCmd(node_rx, RX_RF_ALL, 'rx_gain_mode', 'manual');
            RxGainRF = 2;                  % Rx RF Gain in [1:3]
            RxGainBB = 12;                 % Rx Baseband Gain in [0:31]
            wl_interfaceCmd(node_rx, RX_RF_ALL, 'rx_gains', RxGainRF, RxGainBB);
        end

        % Get parameters from the node
        SAMP_FREQ    = wl_basebandCmd(nodes(1), 'tx_buff_clk_freq');
        Ts           = 1/SAMP_FREQ;

        % We will read the transmitter's maximum I/Q buffer length
        % and assign that value to a temporary variable.
        %
        % NOTE:  We assume that the buffers sizes are the same for all interfaces

        maximum_buffer_len = min(MAX_TX_LEN, wl_basebandCmd(node_tx, TX_RF_VEC, 'tx_buff_max_num_samples'));
        example_mode_string = 'hw';
    else
        % Use same defaults for hardware-dependent params in sim-only version
        maximum_buffer_len = min(MAX_TX_LEN, 2^20);
        SAMP_FREQ           = 40e6;
        example_mode_string = 'sim';
    end

    % Willy : generate two random channels
    H_1 = (randn(1, N_SC) + 1i*randn(1, N_SC));
    H_1 = H_1 ./ abs(H_1);
   % H_1 = H_1/sqrt(2);
    H_2 = (randn(1, N_SC) + 1i*randn(1, N_SC));
    H_2 = H_2 ./ abs(H_2);
   % H_1 = H_1/sqrt(2);

    % Willy : generate w1, w2
    w1 = ones(1,N_SC);
    % w1 * h1 + w2 * h2 = 0
    w2 = -(H_1./ H_2).* w1; 
    w_power = abs(w1).^2 + abs(w2).^2;
    w1 = w1 ./ sqrt(w_power);
    w2 = w2 ./ sqrt(w_power);
    
    %% Define a half-band 2x interpolation filter response
    interp_filt2 = zeros( 1 , 43 );
    interp_filt2([1 3 5 7 9 11 13 15 17 19 21]) = [12 -32 72 -140 252 -422 682 -1086 1778 -3284 10364];
    interp_filt2([23 25 27 29 31 33 35 37 39 41 43]) = interp_filt2(fliplr([1 3 5 7 9 11 13 15 17 19 21]));
    interp_filt2(22) = 16384;
    interp_filt2 = interp_filt2./max(abs(interp_filt2));

    % Define the preamble
    % Note: The STS symbols (short preamble) in the preamble meet the requirements needed by the
    % AGC core at the receiver. Details on the operation of the AGC are
    % available on the wiki: http://warpproject.org/trac/wiki/WARPLab/AGC
    sts_f = zeros( 1 , 64 );
    sts_f( 1 : 27 ) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
    sts_f( 39 : 64 ) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
    sts_t = ifft( sqrt( 13 / 6 ) .* sts_f , 64 );
    sts_t = sts_t( 1 : 16 );

    % Willy : precode with h1 and h2 separately
    sts_t1 = ifft( sqrt( 13 / 6 ) .* sts_f .* H_1 , 64 );
    sts_t1 = sts_t1( 1 : 16 );

    sts_t2 = ifft( sqrt( 13 / 6 ) .* sts_f .* H_2 , 64 );
    sts_t2 = sts_t2( 1 : 16 );

    % LTS for CFO and channel estimation
    lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
    lts_t = ifft( lts_f , 64 );
    
    % Willy : precode with h1 and h2 separately
    lts_t_1 = lts_f .* H_1;
    lts_t_2 = lts_f .* H_2;
    
    lts_t_f1 = ifft( lts_t_1 , 64 );
    lts_t_f2 = ifft( lts_t_2 , 64 );
    
    % Use 30 copies of the 16-sample STS for extra AGC settling margin
    preamble = [repmat(sts_t, 1, 30)  lts_t(33:64) lts_t lts_t];
    
    % Willy : preamble
    preamble1 = [repmat(sts_t1, 1, 30)  lts_t_f1(33:64) lts_t_f1 lts_t_f1];
    preamble2 = [repmat(sts_t2, 1, 30)  lts_t_f2(33:64) lts_t_f2 lts_t_f2];

    % Sanity check variables that affect the number of Tx samples
    num_samps_needed = ceil((TRIGGER_OFFSET_TOL_NS*1e-9) / (1/SAMP_FREQ)) + ...
                        INTERP_RATE*((N_OFDM_SYMS * (N_SC + CP_LEN)) + length(preamble) +  ceil(length(interp_filt2)/2));

    if num_samps_needed > maximum_buffer_len
        fprintf('Too many OFDM symbols for TX_NUM_SAMPS!\n');
        fprintf('Raise MAX_TX_LEN to %d, or \n', num_samps_needed);
        fprintf('Reduce N_OFDM_SYMS to %d\n', floor(((maximum_buffer_len - ceil((TRIGGER_OFFSET_TOL_NS*1e-9) / (1/SAMP_FREQ)))/INTERP_RATE - (length(preamble) +  ceil(length(interp_filt2)/2)))/(N_SC + CP_LEN)));
        return;
    end


    %% Generate a payload of random integers
    % Willy : modulate another sequence 
    tx_data = randi( MOD_ORDER , 1 , N_DATA_SYMS ) - 1;

    % Functions for data -> complex symbol mapping (like qammod, avoids comm toolbox requirement)
    % These anonymous functions implement the modulation mapping from IEEE 802.11-2012 Section 18.3.5.8
    modvec_bpsk   =  1 .* [-1 1];
    modvec_16qam  =  ( 1 / sqrt( 10 ) ) .* [-3 -1 +3 +1];
    modvec_64qam  =  ( 1 / sqrt( 43 ) ) .* [-7 -5 -1 -3 +7 +5 +1 +3];

    mod_fcn_bpsk  = @(x) complex(modvec_bpsk(1+x),0);
    mod_fcn_qpsk  = @(x) complex(modvec_bpsk(1+bitshift(x, -1)), modvec_bpsk(1+mod(x, 2)));
    mod_fcn_16qam = @(x) complex(modvec_16qam(1+bitshift(x, -2)), modvec_16qam(1+mod(x,4)));
    mod_fcn_64qam = @(x) complex(modvec_64qam(1+bitshift(x, -3)), modvec_64qam(1+mod(x,8)));

    % Map the data values on to complex symbols
    switch MOD_ORDER
        case 2         % BPSK
            % Willy : two sequences
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
    tx_syms_mat = reshape( tx_syms , length( SC_IND_DATA ) , N_OFDM_SYMS );

    % Define the pilot tone values as BPSK symbols
    pilots = [ 1 , 1 , -1 , 1 ].';

    % Repeat the pilots across all OFDM symbols
    pilots_mat = repmat( pilots , 1 , N_OFDM_SYMS );


    %% Willy : IFFT : tx1 and tx2

    % Construct the IFFT input matrix
    ifft_in_mat = zeros( N_SC , N_OFDM_SYMS );
    
    % Willy : generate the matrix for precoding and ifft
    ifft_in_mat1 = zeros( N_SC , N_OFDM_SYMS );
    ifft_in_mat2 = zeros( N_SC , N_OFDM_SYMS );

    % Insert the data and pilot values; other subcarriers will remain at 0
    ifft_in_mat( SC_IND_DATA , : ) = tx_syms_mat;
    ifft_in_mat( SC_IND_PILOTS , : ) = pilots_mat;
    
    % Willy : precoding
    ifft_in_mat1( SC_IND_DATA , : ) = tx_syms_mat;
    ifft_in_mat1( SC_IND_PILOTS , : ) = pilots_mat;

    ifft_in_mat2( SC_IND_DATA , : ) = tx_syms_mat;
    ifft_in_mat2( SC_IND_PILOTS , : ) = pilots_mat;

    H_1 = repmat(H_1, N_OFDM_SYMS, 1);
    H_2 = repmat(H_2, N_OFDM_SYMS, 1);
    w1 = repmat(w1, N_OFDM_SYMS, 1);
    w2 = repmat(w2, N_OFDM_SYMS, 1);
    
    % Willy : without nulling
    freq_yn1 = ifft_in_mat1 .* H_1.'; 
    freq_yn2 = ifft_in_mat2 .* H_2.';
    
    % Willy : Perform the IFFT for freq_y1, freq_y2
    tx_payload_matn1 = ifft(freq_yn1, N_SC , 1);
    tx_payload_matn2 = ifft(freq_yn2, N_SC , 1);
    
     % Willy : Insert CP
    if ( CP_LEN > 0 )
        tx_cpn1 = tx_payload_matn1( ( end-CP_LEN + 1 : end ) , : );
        tx_payload_matn1 = [ tx_cpn1 ; tx_payload_matn1 ];
    end

    if ( CP_LEN > 0 )
        tx_cpn2 = tx_payload_matn2( ( end-CP_LEN + 1 : end ) , : );
        tx_payload_matn2 = [ tx_cpn2 ; tx_payload_matn2 ];
    end
    
    % Willy : Reshape
    tx_payload_vecn1 = reshape( tx_payload_matn1 , 1 , numel( tx_payload_matn1 ) );
    tx_payload_vecn2 = reshape( tx_payload_matn2 , 1 , numel( tx_payload_matn2 ) );
    
    % Willy : with nulling
    freq_y1 = freq_yn1 .* w1.';
    freq_y2 = freq_yn2 .* w2.';

    
    % Willy : Perform the IFFT for freq_y1, freq_y2
    tx_payload_mat1 = ifft(freq_y1, N_SC , 1);
    tx_payload_mat2 = ifft(freq_y2, N_SC , 1);
   
    
     % Willy : Insert CP
    if ( CP_LEN > 0 )
        tx_cp1 = tx_payload_mat1( ( end-CP_LEN + 1 : end ) , : );
        tx_payload_mat1 = [ tx_cp1 ; tx_payload_mat1 ];
    end

    if ( CP_LEN > 0 )
        tx_cp2 = tx_payload_mat2( ( end-CP_LEN + 1 : end ) , : );
        tx_payload_mat2 = [ tx_cp2 ; tx_payload_mat2 ];
    end
    
    % Willy : Reshape
    tx_payload_vec1 = reshape( tx_payload_mat1 , 1 , numel( tx_payload_mat1 ) );
    tx_payload_vec2 = reshape( tx_payload_mat2 , 1 , numel( tx_payload_mat2 ) );
    
    % Perform the IFFT
    tx_payload_mat = ifft( ifft_in_mat , N_SC , 1 );

    % Insert the cyclic prefix
    if ( CP_LEN > 0 )
        tx_cp = tx_payload_mat( ( end-CP_LEN + 1 : end ) , : );
        tx_payload_mat = [ tx_cp ; tx_payload_mat ];
    end
    
    % Reshape to a vector
    tx_payload_vec = reshape( tx_payload_mat , 1 , numel( tx_payload_mat ) );

    % Construct the full time-domain OFDM waveform
    tx_vec = [ preamble , tx_payload_vec ];
    
    % Willy : Construct the full time-domain OFDM waveform
    tx_vec1 = [ preamble1 , zeros(1,length(preamble2)) , tx_payload_vecn1 , zeros(1,length(tx_payload_vecn2)) , tx_payload_vecn1 ./ sqrt(2) , tx_payload_vec1 ];
    tx_vec2 = [ zeros(1,length(preamble1)) , preamble2 , zeros(1,length(tx_payload_vecn1)) , tx_payload_vecn2 , tx_payload_vecn2 ./ sqrt(2) , tx_payload_vec2 ];
    
    %tx_vec_air1 = tx_vec1;
    %tx_vec_air2 = tx_vec2;
    scale_factor = max(abs([tx_vec1  tx_vec2]));
    tx_vec1 = TX_SCALE .* tx_vec1 ./ scale_factor;
    tx_vec2 = TX_SCALE .* tx_vec2 ./ scale_factor;
    
    % Pad with zeros for transmission to deal with delay through the
    % interpolation filter
    tx_vec_padded = [ tx_vec , zeros( 1 , ceil( length( interp_filt2 ) / 2 ) ) ];
    
    % Willy : interpolation filter
    tx_vec_padded1 = [ tx_vec1 , zeros( 1 , ceil( length( interp_filt2 ) / 2 ) ) ];
    tx_vec_padded2 = [ tx_vec2 , zeros( 1 , ceil( length( interp_filt2 ) / 2 ) ) ];

    %% Interpolate
    % Zero pad then filter (same as interp or upfirdn without signal processing toolbox)
%{    
    if INTERP_RATE ~= 2
       fprintf( 'Error: INTERP_RATE must equal 2\n' ); 
       return;
    end
    

    tx_vec_2x = zeros( 1 , 2 * numel( tx_vec_padded ) );
    tx_vec_2x( 1 : 2 : end ) = tx_vec_padded;
    tx_vec_air = filter( interp_filt2 , 1 , tx_vec_2x );
    
    % Willy : construct tx_vec_air
    tx_vec_2x1 = zeros( 1 , 2 * numel( tx_vec_padded1 ) );
    tx_vec_2x1( 1 : 2 : end ) = tx_vec_padded1;
    tx_vec_air1 = filter( interp_filt2 , 1 , tx_vec_2x1 );
    
    tx_vec_2x2 = zeros( 1 , 2 * numel( tx_vec_padded2 ) );
    tx_vec_2x2( 1 : 2 : end ) = tx_vec_padded2;
    tx_vec_air2 = filter( interp_filt2 , 1 , tx_vec_2x2 );

    % Scale the Tx vector to +/- 1
    tx_vec_air = TX_SCALE .* tx_vec_air ./ max( abs( tx_vec_air ) );
    % generate in using WARP
    TX_NUM_SAMPS = length( tx_vec_air ); 
    
    % Willy : Scale the Tx vector
    tx_vec_air1 = TX_SCALE .* tx_vec_air1 ./ max( abs( tx_vec_air1 ) );
    tx_vec_air2 = TX_SCALE .* tx_vec_air2 ./ max( abs( tx_vec_air2 ) );

    if ( USE_WARPLAB_TXRX )
        wl_basebandCmd( nodes , 'tx_delay' , 0 );
        % Number of samples to send
        wl_basebandCmd( nodes , 'tx_length' , TX_NUM_SAMPS );
        % Number of samples to receive
        wl_basebandCmd( nodes , 'rx_length' , TX_NUM_SAMPS + ceil( ( TRIGGER_OFFSET_TOL_NS * 1e-9 ) / ( 1 / SAMP_FREQ ) ) );
    end
    
    if USE_WARPLAB_TXRX
        % Write the Tx waveform to the Tx node
        wl_basebandCmd(node_tx, TX_RF_VEC, 'write_IQ', tx_vec_air(:));

        % Enable the Tx and Rx radios
        wl_interfaceCmd(node_tx, TX_RF, 'tx_en');
        wl_interfaceCmd(node_rx, RX_RF, 'rx_en');

        % Enable the Tx and Rx buffers
        wl_basebandCmd(node_tx, TX_RF, 'tx_buff_en');
        wl_basebandCmd(node_rx, RX_RF, 'rx_buff_en');

        % Trigger the Tx/Rx cycle at both nodes
        eth_trig.send();

        % Retrieve the received waveform from the Rx node
        rx_vec_air = wl_basebandCmd(node_rx, RX_RF_VEC, 'read_IQ', 0, TX_NUM_SAMPS + (ceil((TRIGGER_OFFSET_TOL_NS*1e-9) / (1/SAMP_FREQ))));

        rx_vec_air = rx_vec_air(:).';

        % Disable the Tx/Rx radios and buffers
        wl_basebandCmd(node_tx, TX_RF_ALL, 'tx_rx_buff_dis');
        wl_basebandCmd(node_rx, RX_RF_ALL, 'tx_rx_buff_dis');

        wl_interfaceCmd(node_tx, TX_RF_ALL, 'tx_rx_dis');
        wl_interfaceCmd(node_rx, RX_RF_ALL, 'tx_rx_dis');
    else
        % Sim-only mode: Apply wireless degradations here for sim (noise, fading, etc)

        % Perfect (ie. Rx=Tx):
        % rx_vec_air = tx_vec_air;
        % AWGN:
        rx_vec_air = [ tx_vec_air , zeros( 1 , ceil( ( TRIGGER_OFFSET_TOL_NS * 1e-9 ) / ( 1 / SAMP_FREQ ) ) ) ];
        
        % Willy : AWGN
        rx_vec_air1 = [ tx_vec_air1 , zeros( 1 , ceil( ( TRIGGER_OFFSET_TOL_NS * 1e-9 ) / ( 1 / SAMP_FREQ ) ) ) ];
        rx_vec_air2 = [ tx_vec_air2 , zeros( 1 , ceil( ( TRIGGER_OFFSET_TOL_NS * 1e-9 ) / ( 1 / SAMP_FREQ ) ) ) ];
        
        Pn = ( mean( abs( rx_vec_air ) .^ 2 ) ) / SNR;
        %rx_vec_air = rx_vec_air + sqrt( Pn ) * complex( randn( 1 , length( rx_vec_air ) ), randn( 1 , length( rx_vec_air ) ) ) / sqrt( 2 );
        
        % Willy : generate and add the noise 
        %rx_vec_air_mix = rx_vec_air1 + rx_vec_air2;
        rx_vec_air_mix = tx_vec1 + tx_vec2;
        noise_power = mean((abs( tx_vec1 ).^2) + (abs( tx_vec2 ).^2))/ SNR;
        rx_vec_air_mix = rx_vec_air_mix + sqrt( noise_power ) * complex( randn( 1 , length( rx_vec_air_mix ) ), randn( 1 , length( rx_vec_air_mix ) ) ) / sqrt( 2 );
        
        % CFO:
        % rx_vec_air = tx_vec_air .* exp(-1i * 2 * pi * 1e-4 * [0 : length(tx_vec_air) - 1]);
        end
    %}
 % Step 3. Simulated SNR: 5, 10, 15, 20 and 25 dB.
        tx_vec1 = [tx_vec1, zeros(1,ceil((TRIGGER_OFFSET_TOL_NS*1e-9) / (1/SAMP_FREQ)))];
        tx_vec2 = [tx_vec2, zeros(1,ceil((TRIGGER_OFFSET_TOL_NS*1e-9) / (1/SAMP_FREQ)))];
        SNR_dB = SNR_simulated;
        %SNR_dB = 25;
        SNR = 10 ^ ( SNR_dB / 10 );    
% Willy : generate and add the noise 
        %rx_vec_air_mix = rx_vec_air1 + rx_vec_air2;
        rx_vec_air_mix = tx_vec1 + tx_vec2;
        noise_power = mean((abs( tx_vec1 ).^2) + (abs( tx_vec2 ).^2))/ SNR;
        rx_vec_air_mix = rx_vec_air_mix + sqrt( noise_power ) * complex( randn( 1 , length( rx_vec_air_mix ) ), randn( 1 , length( rx_vec_air_mix ) ) ) / sqrt( 2 );
        
    %% Output files
    file = fopen( 'tx_data.bin' , 'w' );
    fwrite( file , tx_data , 'int' );
    fclose( file );
    
    z_real = real( tx_syms_mat );
    z_imag = imag( tx_syms_mat );
    adjacent = [ z_real , z_imag ];
    file = fopen( 'tx_syms_mat.bin' , 'w' );
    fwrite( file , adjacent , 'float' );
    fclose( file );
    
    complex_array1 = zeros( 1 , 2 * length( tx_vec1 ) );
    for j = 1 : length( tx_vec1 )
        complex_array1( 2 * j - 1 ) = real( tx_vec1( j ) );
        complex_array1( 2 * j ) = imag( tx_vec1( j ) );
    end
    
    file = fopen( 'rx_vec_air1.bin' , 'w' );
    fwrite( file , complex_array1 , 'float' );
    fclose( file );
    
    complex_array2 = zeros( 1 , 2 * length( tx_vec2 ) );
    for j = 1 : length( tx_vec2 )
        complex_array2( 2 * j - 1 ) = real( tx_vec2( j ) );
        complex_array2( 2 * j ) = imag( tx_vec2( j ) );
    end
    
    file = fopen( 'rx_vec_air2.bin' , 'w' );
    fwrite( file , complex_array2 , 'float' );
    fclose( file );
    
    
    complex_arraym = zeros( 1 , 2 * length( rx_vec_air_mix ) );
    for j = 1 : length( rx_vec_air_mix )
        complex_arraym( 2 * j - 1 ) = real( rx_vec_air_mix( j ) );
        complex_arraym( 2 * j ) = imag( rx_vec_air_mix( j ) );
    end
    
    file = fopen( 'rx_vec_air_mix.bin' , 'w' );
    fwrite( file , complex_arraym , 'float' );
    fclose( file );
    
end
