function [ SNR_decoded_dB , SNR_decoded_dB1, SNR_decoded_dB2, SNR_decoded_dB3, R_Sig ] = decode( ~ )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % wl_example_siso_ofdm_txrx.m
    % A detailed write-up of this example is available on the wiki:
    % http://warpproject.org/trac/wiki/WARPLab/Examples/OFDM
    %
    % Copyright (c) 2015 Mango Communications - All Rights Reserved
    % Distributed under the WARP License (http://warpproject.org/license)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
%clear;
    %% Rx Site
    % Params:
    % Change into simualtion mode.
    USE_WARPLAB_TXRX        = 0;           % Enable WARPLab-in-the-loop (otherwise sim-only)
    WRITE_PNG_FILES         = 0;           % Enable writing plots to PNG
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

    % Rx processing params
    FFT_OFFSET                    = 0;           % Number of CP samples to use in FFT (on average)
    % Normalized threshold for LTS correlation for packet detection
    LTS_CORR_THRESH               = 0.8;         
    DO_APPLY_CFO_CORRECTION       = 0;           % Enable CFO estimation/correction
    DO_APPLY_PHASE_ERR_CORRECTION = 0;           % Enable Residual CFO estimation/correction
    DO_APPLY_SFO_CORRECTION       = 0;           % Enable SFO estimation/correction
    DECIMATE_RATE                 = INTERP_RATE;

    % WARPLab experiment params
    MAX_TX_LEN              = 2 ^ 20;      % Maximum number of samples to use for this experiment

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

    %% Load files
    fid_1 = fopen( 'tx_data.bin','rb' );
    tx_data = fread( fid_1, 'float' );
    tx_data = tx_data.';
    fclose( fid_1 );
    
    file = fopen( 'tx_syms_mat.bin' );
    same_real = fread( file , [ 48 , 50 ] , 'float' );
    same_imag = fread( file , [ 48 , 50 ] , 'float' );
    fclose( file );
    tx_syms_mat = complex( same_real , same_imag );
    %tx_syms_mat = read_complex_binary( 'tx_syms_mat.bin' ).';

    rx_vec_air_mix = read_complex_binary( 'rx_vec_air_mix.bin' );
    R_Sig = rx_vec_air_mix;
    %{
    figure();
    plot(abs(rx_vec_air_mix));
    title("Y Vec Air Mix");
    %}
 % Define the pilot tone values as BPSK symbols
    pilots = [ 1 , 1 , -1 , 1 ].';

    % Repeat the pilots across all OFDM symbols
    pilots_mat = repmat( pilots , 1 , N_OFDM_SYMS );

    lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
    %% Decimate
    if INTERP_RATE ~= 2
       fprintf( 'Error: INTERP_RATE must equal 2\n' ); 
       return;
    end

    %raw_rx_dec = filter( interp_filt2 , 1 , rx_vec_air_mix );
    %raw_rx_dec = raw_rx_dec( 1 : 2 : end );
    raw_rx_dec = rx_vec_air_mix;

    %% Correlate for LTS
%{
    % Complex cross correlation of Rx waveform with time-domain LTS
    lts_corr = abs( conv( conj( fliplr( lts_t ) ) , sign( raw_rx_dec ) ) );

    % Skip early and late samples - avoids occasional false positives from pre-AGC samples
    lts_corr = lts_corr( 32 : end - 32 );

    % Find all correlation peaks
    lts_peaks = find( lts_corr( 1 : 800 ) > LTS_CORR_THRESH * max( lts_corr ) );

    % Select best candidate correlation peak as LTS-payload boundary
    [ LTS1 , LTS2 ] = meshgrid( lts_peaks , lts_peaks );
    [ lts_second_peak_index , y ] = find( LTS2 - LTS1 == length( lts_t ) );

    % Stop if no valid correlation peak was found
    if isempty( lts_second_peak_index )
        fprintf( 'No LTS Correlation Peaks Found!\n' );
        return;
    end

    % Set the sample indices of the payload symbols and preamble
    % The "+32" corresponds to the 32-sample cyclic prefix on the preamble LTS
    % The "-160" corresponds to the length of the preamble LTS (2.5 copies of 64-sample LTS)
    payload_ind = lts_peaks( max( lts_second_peak_index ) ) + 32;
   %}
    payload_ind = 640 * 2 + 1;
    lts_ind = payload_ind - 640 - 160;
    lts_ind2 = lts_ind + 640;
    
    % Willy : lts_ind: H1
    if DO_APPLY_CFO_CORRECTION
        %Extract LTS (not yet CFO corrected)
        rx_lts = raw_rx_dec( lts_ind : lts_ind + 159 );
        rx_lts1 = rx_lts( -64 + -FFT_OFFSET + [ 97 : 160 ] );
        rx_lts2 = rx_lts( -FFT_OFFSET + [ 97 : 160 ]);

        %Calculate coarse CFO est
        rx_cfo_est_lts = mean( unwrap( angle( rx_lts2 .* conj( rx_lts1 ) ) ) );
        rx_cfo_est_lts = rx_cfo_est_lts / ( 2 * pi * 64 );
    else
        rx_cfo_est_lts = 0;
    end

    % Apply CFO correction to raw Rx waveform
    rx_cfo_corr_t = exp( -1i * 2 * pi * rx_cfo_est_lts * [ 0 : length( raw_rx_dec ) - 1 ] );
    rx_dec_cfo_corr1 = raw_rx_dec .* rx_cfo_corr_t;

    % Re-extract LTS for channel estimate
    rx_lts = raw_rx_dec( lts_ind : lts_ind + 159 );
    rx_lts1 = rx_lts( -64 + -FFT_OFFSET + [ 97 : 160 ] );
    rx_lts2 = rx_lts( -FFT_OFFSET + [ 97 : 160 ] );

    rx_lts1_f = fft( rx_lts1 );
    rx_lts2_f = fft( rx_lts2 );

    % Calculate channel estimate from average of 2 training symbols
    rx_H1_est = lts_f .* ( rx_lts1_f + rx_lts2_f ) / 2;

    % Willy : lts_ind2: H2
    if DO_APPLY_CFO_CORRECTION
        %Extract LTS (not yet CFO corrected)
        rx_lts = raw_rx_dec( lts_ind2 : lts_ind2 + 159 );
        rx_lts1 = rx_lts( -64 + -FFT_OFFSET + [ 97 : 160 ] );
        rx_lts2 = rx_lts( -FFT_OFFSET + [ 97 : 160 ]);

        %Calculate coarse CFO est
        rx_cfo_est_lts = mean( unwrap( angle( rx_lts2 .* conj( rx_lts1 ) ) ) );
        rx_cfo_est_lts = rx_cfo_est_lts / ( 2 * pi * 64 );
    else
        rx_cfo_est_lts = 0;
    end

    % Apply CFO correction to raw Rx waveform
    rx_cfo_corr_t = exp( -1i * 2 * pi * rx_cfo_est_lts * [ 0 : length( raw_rx_dec ) - 1 ] );
    rx_dec_cfo_corr2 = raw_rx_dec .* rx_cfo_corr_t;

    % Re-extract LTS for channel estimate
    rx_lts = raw_rx_dec( lts_ind2 : lts_ind2 + 159 );
    rx_lts11 = rx_lts( -64 + -FFT_OFFSET + [ 97 : 160 ] );
    rx_lts22 = rx_lts( -FFT_OFFSET + [ 97 : 160 ] );

    rx_lts1_f1 = fft( rx_lts11 );
    rx_lts2_f2 = fft( rx_lts22 );

    % Calculate channel estimate from average of 2 training symbols
    rx_H2_est = lts_f .* ( rx_lts1_f1 + rx_lts2_f2 ) / 2;
    %% Rx payload processing

    % Extract the payload samples (integral number of OFDM symbols following preamble)
    %payload_vec = rx_dec_cfo_corr1(payload_ind : payload_ind+N_OFDM_SYMS*(N_SC+CP_LEN)-1);
    payload_vec = rx_dec_cfo_corr1(payload_ind : payload_ind+4*N_OFDM_SYMS*(N_SC+CP_LEN)-1);
    payload_mat = reshape(payload_vec, (N_SC+CP_LEN), 4*N_OFDM_SYMS);

    % Remove the cyclic prefix, keeping FFT_OFFSET samples of CP (on average)
    payload_mat_noCP = payload_mat(CP_LEN-FFT_OFFSET+[1:N_SC], :);
   
    % Take the FFT
    syms_f_mat = fft(payload_mat_noCP, N_SC, 1);

    % Equalize (zero-forcing, just divide by complex channel estimates)
    %syms_eq_mat = syms_f_mat ./ repmat(rx_H_est.', 1, N_OFDM_SYMS);
    syms_eq_mat = syms_f_mat(:, 1: N_OFDM_SYMS) ./ repmat(rx_H1_est.', 1, N_OFDM_SYMS);
    syms_eq_mat1 = syms_f_mat(:, N_OFDM_SYMS + 1: 2*N_OFDM_SYMS) ./ repmat(rx_H2_est.', 1, N_OFDM_SYMS);
    syms_eq_mat2 = syms_f_mat(:, 2*N_OFDM_SYMS + 1: 3*N_OFDM_SYMS) ./ repmat(rx_H2_est.'+rx_H1_est.', 1, N_OFDM_SYMS);
    syms_eq_mat3 = syms_f_mat(:, 3*N_OFDM_SYMS + 1: 4*N_OFDM_SYMS) ./ repmat(rx_H2_est.'+rx_H1_est.', 1, N_OFDM_SYMS);

    if DO_APPLY_SFO_CORRECTION
        % SFO manifests as a frequency-dependent phase whose slope increases
        % over time as the Tx and Rx sample streams drift apart from one
        % another. To correct for this effect, we calculate this phase slope at
        % each OFDM symbol using the pilot tones and use this slope to
        % interpolate a phase correction for each data-bearing subcarrier.

        % Extract the pilot tones and "equalize" them by their nominal Tx values
        pilots_f_mat = syms_eq_mat( SC_IND_PILOTS , : );
        pilots_f_mat2 = syms_eq_mat1( SC_IND_PILOTS , : );
        pilots_f_mat3 = syms_eq_mat2( SC_IND_PILOTS , : );
        pilots_f_mat4 = syms_eq_mat3( SC_IND_PILOTS , : );
        
        pilots_f_mat_comp = pilots_f_mat .* pilots_mat;
        pilots_f_mat_comp1 = pilots_f_mat1 .* pilots_mat;
        pilots_f_mat_comp2 = pilots_f_mat2 .* pilots_mat;
        pilots_f_mat_comp3 = pilots_f_mat3 .* pilots_mat;

        % Calculate the phases of every Rx pilot tone
        pilot_phases = unwrap(angle(fftshift(pilots_f_mat_comp,1)), [], 1);
        pilot_phases1 = unwrap(angle(fftshift(pilots_f_mat_comp1,1)), [], 1);
        pilot_phases2 = unwrap(angle(fftshift(pilots_f_mat_comp2,1)), [], 1);
        pilot_phases3 = unwrap(angle(fftshift(pilots_f_mat_comp3,1)), [], 1);

        % Calculate slope of pilot tone phases vs frequency in each OFDM symbol
        pilot_spacing_mat = repmat(mod(diff(fftshift(SC_IND_PILOTS)),64).', 1, N_OFDM_SYMS);                        
        pilot_slope_mat = mean(diff(pilot_phases) ./ pilot_spacing_mat);

        % Calculate the SFO correction phases for each OFDM symbol
        
        pilot_phase_sfo_corr = fftshift((-32:31).' * pilot_slope_mat, 1);
        pilot_phase_corr = exp(-1i*(pilot_phase_sfo_corr));

        % Apply the pilot phase correction per symbol
        syms_eq_mat = syms_eq_mat .* pilot_phase_corr;
    else
        % Define an empty SFO correction matrix (used by plotting code below)
        pilot_phase_sfo_corr = zeros(N_SC, N_OFDM_SYMS);
    end


    if DO_APPLY_PHASE_ERR_CORRECTION
        % Extract the pilots and calculate per-symbol phase error
        pilots_f_mat = syms_eq_mat(SC_IND_PILOTS, :);
        pilots_f_mat_comp = pilots_f_mat.*pilots_mat;
        pilot_phase_err = angle(mean(pilots_f_mat_comp));
    else
        % Define an empty phase correction vector (used by plotting code below)
        pilot_phase_err = zeros(1, N_OFDM_SYMS);
    end
    pilot_phase_err_corr = repmat( pilot_phase_err , N_SC , 1 );
    pilot_phase_corr = exp( -1i * ( pilot_phase_err_corr ) );

    % Apply the pilot phase correction per symbol
    syms_eq_pc_mat = syms_eq_mat .* pilot_phase_corr;
    syms_eq_pc_mat1 = syms_eq_mat1 .* pilot_phase_corr;
    syms_eq_pc_mat2 = syms_eq_mat2 .* pilot_phase_corr;
    syms_eq_pc_mat3 = syms_eq_mat3 .* pilot_phase_corr;
    
    payload_syms_mat = syms_eq_pc_mat( SC_IND_DATA , : );
    payload_syms_mat1 = syms_eq_pc_mat1( SC_IND_DATA , : );
    payload_syms_mat2 = syms_eq_pc_mat2( SC_IND_DATA , : );
    payload_syms_mat3 = syms_eq_pc_mat3( SC_IND_DATA , : );

    %% Demodulate
    rx_syms = reshape( payload_syms_mat , 1 , N_DATA_SYMS );
    rx_syms1 = reshape( payload_syms_mat1 , 1 , N_DATA_SYMS );
    rx_syms2 = reshape( payload_syms_mat2 , 1 , N_DATA_SYMS );
    rx_syms3 = reshape( payload_syms_mat3 , 1 , N_DATA_SYMS );

    demod_fcn_bpsk = @(x) double(real(x)>0);
    demod_fcn_qpsk = @(x) double(2*(real(x)>0) + 1*(imag(x)>0));
    demod_fcn_16qam = @(x) (8*(real(x)>0)) + (4*(abs(real(x))<0.6325)) + (2*(imag(x)>0)) + (1*(abs(imag(x))<0.6325));
    demod_fcn_64qam = @(x) (32*(real(x)>0)) + (16*(abs(real(x))<0.6172)) + (8*((abs(real(x))<(0.9258))&&((abs(real(x))>(0.3086))))) + (4*(imag(x)>0)) + (2*(abs(imag(x))<0.6172)) + (1*((abs(imag(x))<(0.9258))&&((abs(imag(x))>(0.3086)))));

    % Willy : Calculate for SNR
    N0 = mean( abs( payload_syms_mat - tx_syms_mat ) .^ 2 );
    SNR_decoded = mean( abs( tx_syms_mat ) .^ 2 ) / N0;
    SNR_decoded_dB = 10 * log10( SNR_decoded );
    
    N0 = mean( abs( payload_syms_mat1 - tx_syms_mat ) .^ 2 );
    SNR_decoded = mean( abs( tx_syms_mat ) .^ 2 ) / N0;
    SNR_decoded_dB1 = 10 * log10( SNR_decoded );
    
    N0 = mean( abs( payload_syms_mat2 - tx_syms_mat ) .^ 2 );
    SNR_decoded = mean( abs( tx_syms_mat ) .^ 2 ) / N0;
    SNR_decoded_dB2 = 10 * log10( SNR_decoded );
    
    N0 = mean( abs( payload_syms_mat3 - tx_syms_mat ) .^ 2 );
    SNR_decoded = mean( abs( tx_syms_mat ) .^ 2 ) / N0;
    SNR_decoded_dB3 = 10 * log10( SNR_decoded );

    switch(MOD_ORDER)
        case 2         % BPSK
            rx_data = arrayfun(demod_fcn_bpsk, rx_syms);
        case 4         % QPSK
            rx_data = arrayfun(demod_fcn_qpsk, rx_syms);
        case 16        % 16-QAM
            rx_data = arrayfun(demod_fcn_16qam, rx_syms);
        case 64        % 64-QAM
            rx_data = arrayfun(demod_fcn_64qam, rx_syms);
    end

    %% Plot Results
    cf = 0;
%{
    % Tx signal
    cf = cf + 1;
    figure(cf); 
    clf;

    subplot( 2 , 1 , 1 );
    plot( real(tx_vec_air), 'b');
    axis( [ 0 , length(tx_vec_air) , -TX_SCALE TX_SCALE ] )
    grid on;
    title('Tx Waveform (I)' );

    subplot( 2 , 1 , 2 );
    plot( imag( tx_vec_air ) , 'r' );
    axis( [ 0 , length( tx_vec_air ) , -TX_SCALE TX_SCALE ] );
    grid on;
    title( 'Tx Waveform (Q)' );

    if WRITE_PNG_FILES
        print( gcf, sprintf( 'wl_ofdm_plots_%s_txIQ' , example_mode_string ) , '-dpng' , '-r96' , '-painters' );
    end

    % Rx signal
    cf = cf + 1;
    figure( cf); 
    clf;
    subplot( 2 , 1 , 1 );
    plot( real( rx_vec_air ) , 'b' );
    axis( [ 0 , length( rx_vec_air ) , -TX_SCALE , TX_SCALE ] );
    grid on;
    title( 'Rx Waveform (I)' );

    subplot( 2 , 1 , 2 );
    plot( imag( rx_vec_air ) , 'r' );
    axis( [ 0 , length( rx_vec_air ) , -TX_SCALE , TX_SCALE ] );
    grid on;
    title( 'Rx Waveform (Q)' );

    if WRITE_PNG_FILES
        print( gcf , sprintf( 'wl_ofdm_plots_%s_rxIQ' , example_mode_string ) , '-dpng' , '-r96' , '-painters' );
    end

    % Rx LTS correlation
    cf = cf + 1;
    figure( cf ); 
    clf;

    lts_to_plot = lts_corr;
    plot( lts_to_plot , '.-b' , 'LineWidth' , 1 );
    hold on;
    grid on;
    line( [ 1 , length( lts_to_plot ) ] , LTS_CORR_THRESH * max( lts_to_plot ) * [ 1 , 1 ] , 'LineStyle' , '--' , 'Color' , 'r' , 'LineWidth' , 2 );
    title( 'LTS Correlation and Threshold' );
    xlabel( 'Sample Index' );
    myAxis = axis();
    axis( [ 1 , 1000 , myAxis( 3 ) , myAxis( 4 ) ] );

    if WRITE_PNG_FILES
        print( gcf , sprintf( 'wl_ofdm_plots_%s_ltsCorr' , example_mode_string ), '-dpng' , '-r96' , '-painters' );
    end

    % Channel Estimates
    cf = cf + 1;

    rx_H_est_plot = repmat( complex( NaN , NaN ) , 1 , length( rx_H_est ) );
    rx_H_est_plot( SC_IND_DATA ) = rx_H_est( SC_IND_DATA );
    rx_H_est_plot( SC_IND_PILOTS ) = rx_H_est( SC_IND_PILOTS );

    x = ( 20 / N_SC ) * ( -( N_SC / 2 ) : ( N_SC / 2 - 1 ) );

    figure( cf ); 
    clf;
    subplot( 2 , 1 , 1 );
    stairs( x - ( 20 / ( 2 * N_SC ) ) , fftshift( real( rx_H_est_plot ) ) , 'b' , 'LineWidth' , 2 );
    hold on;
    stairs( x - ( 20 / ( 2 * N_SC ) ) , fftshift( imag( rx_H_est_plot ) ) , 'r' , 'LineWidth' , 2 );
    hold off;
    axis( [ min( x ) , max( x ) , -1.1 * max( abs( rx_H_est_plot ) ) , 1.1 * max( abs( rx_H_est_plot ) ) ] );
    grid on;
    title( 'Channel Estimates (I and Q)' );

    subplot( 2 , 1 , 2 );
    bh = bar( x, fftshift( abs( rx_H_est_plot ) ) , 1 , 'LineWidth' , 1 );
    shading flat;
    set( bh , 'FaceColor' , [ 0 , 0 , 1 ] );
    axis( [ min( x ) , max( x ) , 0 , 1.1 * max( abs( rx_H_est_plot ) ) ] );
    grid on;
    title( 'Channel Estimates (Magnitude)' );
    xlabel( 'Baseband Frequency (MHz)' );

    if WRITE_PNG_FILES
        print( gcf , sprintf( 'wl_ofdm_plots_%s_chanEst' , example_mode_string ) , '-dpng', '-r96', '-painters' );
    end

    %% Pilot phase error estimate
    cf = cf + 1;
    figure( cf ); 
    clf;
    subplot( 2 , 1 , 1 )
    plot( pilot_phase_err , 'b', 'LineWidth' , 2 );
    title( 'Phase Error Estimates' );
    xlabel( 'OFDM Symbol Index' );
    ylabel( 'Radians' );
    axis( [ 1 , N_OFDM_SYMS , -3.2 , 3.2 ] );
    grid on

    h = colorbar;
    set( h , 'Visible' , 'off' );

    subplot( 2 , 1 , 2 );
    imagesc( 1 : N_OFDM_SYMS , ( SC_IND_DATA - N_SC / 2 ) , fftshift( pilot_phase_sfo_corr , 1 ) );
    xlabel( 'OFDM Symbol Index' );
    ylabel( 'Subcarrier Index' );
    title( 'Phase Correction for SFO' );
    colorbar;
    myAxis = caxis();
    if myAxis( 2 ) - myAxis( 1 ) < pi
       caxis( [ -pi / 2 , pi / 2 ] );
    end

    if WRITE_PNG_FILES
        print( gcf , sprintf('wl_ofdm_plots_%s_phaseError', example_mode_string), '-dpng', '-r96', '-painters' );
    end

    %% Symbol constellation
    cf = cf + 1;
    figure( cf ); 
    clf;

    plot( payload_syms_mat( : ) , 'ro' , 'MarkerSize' , 1 );
    axis square; 
    axis( 1.5 * [ -1 , 1 , -1 , 1 ] );
    grid on;
    hold on;

    plot( tx_syms_mat( : ) , 'bo' );
    title( 'Tx and Rx Constellations' );
    legend( 'Rx' , 'Tx' , 'Location' , 'EastOutside' );

    if WRITE_PNG_FILES
        print( gcf , sprintf( 'wl_ofdm_plots_%s_constellations' , example_mode_string ) , '-dpng' , '-r96' , '-painters' );
    end
%}
    % EVM & SNR

    %% Calculate Rx stats
%{
    sym_errs = sum( tx_data.' ~= rx_data );
    bit_errs = length( find( dec2bin( bitxor( tx_data.' , rx_data ) , 8 ) == '1' ) );
    rx_evm   = sqrt( sum( ( real( rx_syms ) - real( tx_syms ) ) .^ 2 + ( imag( rx_syms ) - imag( tx_syms ) ) .^ 2 ) / ( length( SC_IND_DATA ) * N_OFDM_SYMS ) );

    % Step 4 - Calculate BER.
    BER_decoded = bit_errs / ( length( SC_IND_DATA ) * N_OFDM_SYMS );

    fprintf( '\nResults:\n' );
    fprintf( 'Num Bytes:   %d\n' , N_DATA_SYMS * log2( MOD_ORDER ) / 8);
    fprintf( 'Sym Errors:  %d (of %d total symbols)\n' , sym_errs , N_DATA_SYMS);
    fprintf( 'Bit Errors:  %d (of %d total bits)\n' , bit_errs , N_DATA_SYMS * log2( MOD_ORDER ) );

    cfo_est_lts = rx_cfo_est_lts * ( SAMP_FREQ / INTERP_RATE );
    cfo_est_phaseErr = mean( diff( unwrap( pilot_phase_err ) ) ) / ( 4e-6 * 2 * pi );
    cfo_total_ppm = ( ( cfo_est_lts + cfo_est_phaseErr ) / ( ( 2.412 + ( .005 * ( CHANNEL - 1 ) ) ) * 1e9 ) ) * 1e6;

    fprintf( 'CFO Est:     %3.2f kHz (%3.2f ppm)\n' , ( cfo_est_lts + cfo_est_phaseErr ) * 1e-3, cfo_total_ppm );
    fprintf( '     LTS CFO Est:                  %3.2f kHz\n' , cfo_est_lts * 1e-3 );
    fprintf( '     Phase Error Residual CFO Est: %3.2f kHz\n' , cfo_est_phaseErr * 1e-3 );

    if DO_APPLY_SFO_CORRECTION
        drift_sec = pilot_slope_mat / ( 2 * pi * 312500 );
        sfo_est_ppm =  1e6 * mean( ( diff( drift_sec ) / 4e-6 ) );
        sfo_est = sfo_est_ppm * 20;
        fprintf( 'SFO Est:     %3.2f Hz (%3.2f ppm)\n' , sfo_est , sfo_est_ppm );
    end
   %}
end