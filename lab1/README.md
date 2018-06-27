# Wireless Communication Systems Lab1 OFDM over USRP

## File structure
* File `single_tx.cpp`          - the TX running on the USRP
* File `single_tx.h`            - the `single_tx.cpp` including file (library)
* File `single_rx.cpp`          - the RX running on the USRP
* File `single_rx.h`            - the `single_rx.cpp` including file (library)
* File `SFO_figure.m`           - script for calling function decode, using to plot the with SFO and without SFO
* File `sigal_gen.m`            - Script of TX, simulating to generate the signal
* File `decode.m`               - function of RX, simulating to decode the signal and plot the figure
* File `read_complex_binary.m`  - function to read complex number
* File `Report_G2.pdf`          - Report of the lab (Group 2) , including the TASK-5 (figure of result, Calculating the SNR and BER) and discussion / ovservation
 
## Usage
1. put the `single_tx.cpp`, `single_tx.h`, `single_rx.cpp`, `single_rx.h` files on the Server
2. modify the `CMakeLists.txt` and use make to compile it
3. use script file `sigal_gen.m` to generate the `tx_vec.bin`
4. run the `single_tx` and `single_rx`
5. download the `.bin` file 
6. use script file `SFO_figure.m` to plot figure, it will run the `decode` function
 

