%% data_Calibration_and_EQ.m
%
% This script finds the required EQ/Calibration transfer functions for the
% measurements
%
% Last updated by Michael Klunk 02/16/2021

%% Initializing the Script

clear
close all
clc

%% Loading in the Calibration/EQ Data

% Individual calibration data

    % PHASE MATCH 1
    load('/Users/michaelklunk/Documents/Msc_Dissertation_Files/Lab Day 2/Mini Mic Phase Match 10-Feb-2022 15:56:22.mat');
    phase_Match_1{1} = data(1).IR;
    phase_Match_1{2} = data(2).IR;
    clear app data specs
    
    % PHASE MATCH 2
    load('/Users/michaelklunk/Documents/Msc_Dissertation_Files/Lab Day 2/Mini Mic Phase Match TWO 10-Feb-2022 15:56:22.mat');
    phase_Match_2{1} = data(1).IR;
    phase_Match_2{2} = data(2).IR;
    clear app data specs
    
% Mini microphone to H.A.T.S. microphone calibration data

    % MINI TO H.A.T.S. CALIBRATION 
    load('/Users/michaelklunk/Documents/Msc_Dissertation_Files/Lab Day 2/Mini to HATS Calibration 10-Feb-2022 16:17:56.mat');
    mini_HATS{1} = data(1).IR;
    mini_HATS{2} = data(2).IR;
    mini_HATS{3} = data(3).IR;
    mini_HATS{4} = data(4).IR;
    clear app data specs
    
%% Organizing the Individual Impulse Responses

% IRs for phase matching calibration

    % LEFT MINI MICROPHONE (PM -> Phase Match)
    IR_PM_Mini_L(:,1) = phase_Match_1{1}(:,1);
    IR_PM_Mini_L(:,2) = phase_Match_1{2}(:,1);
    IR_PM_Mini_L(:,3) = phase_Match_2{1}(:,1);
    IR_PM_Mini_L(:,4) = phase_Match_2{2}(:,1);
    
    % RIGHT MINI MICROPHONE (PM -> Phase Match)
    IR_PM_Mini_R(:,1) = phase_Match_1{1}(:,2);
    IR_PM_Mini_R(:,2) = phase_Match_1{2}(:,2);
    IR_PM_Mini_R(:,3) = phase_Match_2{1}(:,2);
    IR_PM_Mini_R(:,4) = phase_Match_2{2}(:,2);
    
    % LEFT H.A.T.S (PM -> Phase Match)
    IR_PM_HATS_L(:,1) = phase_Match_1{1}(:,3);
    IR_PM_HATS_L(:,2) = phase_Match_1{2}(:,3);
    IR_PM_HATS_L(:,3) = phase_Match_2{1}(:,3);
    IR_PM_HATS_L(:,4) = phase_Match_2{2}(:,3);
    
    % RIGHT H.A.T.S (PM -> Phase Match)
    IR_PM_HATS_R(:,1) = phase_Match_1{1}(:,4);
    IR_PM_HATS_R(:,2) = phase_Match_1{2}(:,4);
    IR_PM_HATS_R(:,3) = phase_Match_2{1}(:,4);
    IR_PM_HATS_R(:,4) = phase_Match_2{2}(:,4);
    
% IRs for microphone relative transfer function calibration

    % LEFT MINI MICROPHONE (MRTF -> Microphone Relative Transfer Function)
    IR_MRTF_Mini_L(:,1) = mini_HATS{1}(:,1);
    IR_MRTF_Mini_L(:,2) = mini_HATS{2}(:,1);
    IR_MRTF_Mini_L(:,3) = mini_HATS{3}(:,1);
    IR_MRTF_Mini_L(:,4) = mini_HATS{4}(:,1);
    
    % RIGHT MINI MICROPHONE (MRTF -> Microphone Relative Transfer Function)
    IR_MRTF_Mini_R(:,1) = mini_HATS{1}(:,2);
    IR_MRTF_Mini_R(:,2) = mini_HATS{2}(:,2);
    IR_MRTF_Mini_R(:,3) = mini_HATS{3}(:,2);
    IR_MRTF_Mini_R(:,4) = mini_HATS{4}(:,2);
    
    % LEFT H.A.T.S (MRTF -> Microphone Relative Transfer Function)
    IR_MRTF_HATS_L(:,1) = mini_HATS{1}(:,3);
    IR_MRTF_HATS_L(:,2) = mini_HATS{2}(:,3);
    IR_MRTF_HATS_L(:,3) = mini_HATS{3}(:,3);
    IR_MRTF_HATS_L(:,4) = mini_HATS{4}(:,3);
    
    % RIGHT H.A.T.S (MRTF -> Microphone Relative Transfer Function)
    IR_MRTF_HATS_R(:,1) = mini_HATS{1}(:,4);
    IR_MRTF_HATS_R(:,2) = mini_HATS{2}(:,4);
    IR_MRTF_HATS_R(:,3) = mini_HATS{3}(:,4);
    IR_MRTF_HATS_R(:,4) = mini_HATS{4}(:,4);

%% Plotting the IRs for the Phase Matching

% Mini Microphones

    % IR MINI LEFT PLOTS
    figure
    subplot(2,2,1)
    plot(IR_PM_Mini_L(:,1))
    hold on
    plot(IR_PM_Mini_L(:,2))
    plot(IR_PM_Mini_L(:,3))
    plot(IR_PM_Mini_L(:,4))
    xlabel('Samples')
    ylabel('Amplitude')
    title('IRs for the Left Mini Microphone (Phase Matching Calibration)')
    legend('First Measurement', 'Second Measurement', 'Third Measurement',...
        'Forth Measurement')
    grid on, grid minor
    
    % IR MINI RIGHT PLOTS
    subplot(2,2,2)
    plot(IR_PM_Mini_R(:,1))
    hold on
    plot(IR_PM_Mini_R(:,2))
    plot(IR_PM_Mini_R(:,3))
    plot(IR_PM_Mini_R(:,4))
    xlabel('Samples')
    ylabel('Amplitude')
    title('IRs for the Right Mini Microphone (Phase Matching Calibration)')
    legend('First Measurement', 'Second Measurement', 'Third Measurement',...
        'Forth Measurement')
    grid on, grid minor
    
% HATS Microphones

    % IR HATS LEFT PLOTS
    subplot(2,2,3)
    plot(IR_PM_HATS_L(:,1))
    hold on
    plot(IR_PM_HATS_L(:,2))
    plot(IR_PM_HATS_L(:,3))
    plot(IR_PM_HATS_L(:,4))
    xlabel('Samples')
    ylabel('Amplitude')
    title('IRs for the Left H.A.T.S. Microphone (Phase Matching Calibration)')
    legend('First Measurement', 'Second Measurement', 'Third Measurement',...
        'Forth Measurement')
    grid on, grid minor
    
    % IR HATS RIGHT PLOTS
    subplot(2,2,4)
    plot(IR_PM_HATS_R(:,1))
    hold on
    plot(IR_PM_HATS_R(:,2))
    plot(IR_PM_HATS_R(:,3))
    plot(IR_PM_HATS_R(:,4))
    xlabel('Samples')
    ylabel('Amplitude')
    title('IRs for the Right H.A.T.S. Microphone (Phase Matching Calibration)')
    legend('First Measurement', 'Second Measurement', 'Third Measurement',...
        'Forth Measurement')
    grid on, grid minor
    
%% Plotting the IRs for the Mini to H.A.T.S. Calibration

% Mini Microphones

    % IR MINI LEFT PLOTS
    figure
    subplot(2,2,1)
    plot(IR_MRTF_Mini_L(:,1))
    hold on
    plot(IR_MRTF_Mini_L(:,2))
    plot(IR_MRTF_Mini_L(:,3))
    plot(IR_MRTF_Mini_L(:,4))
    xlabel('Samples')
    ylabel('Amplitude')
    title('IRs for the Left Mini Microphone (Relative Transfer Function Calibration)')
    legend('First Measurement', 'Second Measurement', 'Third Measurement',...
        'Forth Measurement')
    grid on, grid minor
    
    % IR MINI RIGHT PLOTS
    subplot(2,2,2)
    plot(IR_MRTF_Mini_R(:,1))
    hold on
    plot(IR_MRTF_Mini_R(:,2))
    plot(IR_MRTF_Mini_R(:,3))
    plot(IR_MRTF_Mini_R(:,4))
    xlabel('Samples')
    ylabel('Amplitude')
    title('IRs for the Right Mini Microphone (Relative Transfer Function Calibration)')
    legend('First Measurement', 'Second Measurement', 'Third Measurement',...
        'Forth Measurement')
    grid on, grid minor
    
% HATS Microphones

    % IR HATS LEFT PLOTS
    subplot(2,2,3)
    plot(IR_MRTF_HATS_L(:,1))
    hold on
    plot(IR_MRTF_HATS_L(:,2))
    plot(IR_MRTF_HATS_L(:,3))
    plot(IR_MRTF_HATS_L(:,4))
    xlabel('Samples')
    ylabel('Amplitude')
    title('IRs for the Left H.A.T.S. Microphone (Relative Transfer Function Calibration)')
    legend('First Measurement', 'Second Measurement', 'Third Measurement',...
        'Forth Measurement')
    grid on, grid minor
    
    % IR HATS RIGHT PLOTS
    subplot(2,2,4)
    plot(IR_MRTF_HATS_R(:,1))
    hold on
    plot(IR_MRTF_HATS_R(:,2))
    plot(IR_MRTF_HATS_R(:,3))
    plot(IR_MRTF_HATS_R(:,4))
    xlabel('Samples')
    ylabel('Amplitude')
    title('IRs for the Right H.A.T.S. Microphone (Relative Transfer Function Calibration)')
    legend('First Measurement', 'Second Measurement', 'Third Measurement',...
        'Forth Measurement')
    grid on, grid minor

%% Finding the Calibration Transfer Functions (Phase Matching)
    
% Averaging the respective impulse response sets
    
    % LEFT MINI MICROPHONE IMPULSES
    IR_PM_Mini_L_AVG = mean(IR_PM_Mini_L,2);

    % RIGHT MINI MICROPHONE IMPULSES
    IR_PM_Mini_R_AVG = mean(IR_PM_Mini_R,2);

    % RIGHT H.A.T.S. MICROPHONE IMPULSES
    IR_PM_HATS_L_AVG = mean(IR_PM_HATS_L,2);
    
    % RIGHT H.A.T.S. MICROPHONE IMPULSES
    IR_PM_HATS_R_AVG = mean(IR_PM_HATS_R,2);
    
% Finding the DFTs of the different calibration transfer functions

    % FFT Parameters
    Fs = 48000;                 % Sampling frequency (Hz)
    N = 2^8;                    % Number of FFT Points
    dt = 1/Fs;                  % Delta t (s)
    T = dt*N;                   % Sampling Period (s)
    df = 1/T;                   % Delta f (Hz)
    f = (0:((N/2)-0))./N*Fs;    % Frequency vector (Hz)

    % LEFT MINI MICROPHONE TRANSFER FUNCTION 
    TF_PM_Mini_L = fft(IR_PM_Mini_L_AVG, N);
    
    % RIGHT MINI MICROPHONE TRANSFER FUNCTION 
    TF_PM_Mini_R = fft(IR_PM_Mini_R_AVG, N);
    
    % LEFT H.A.T.S. MICROPHONE TRANSFER FUNCTION 
    TF_PM_HATS_L = fft(IR_PM_HATS_L_AVG, N);
    
    % RIGHT H.A.T.S. MICROPHONE TRANSFER FUNCTION 
    TF_PM_HATS_R = fft(IR_PM_HATS_R_AVG, N);
    
% Calibration transfer functions

    % (mini_Mic_L./mini_Mic_R)
    mini_PM_Calibration = TF_PM_Mini_L./TF_PM_Mini_R;
    
    % (HATS_L./HATS_R)
    HATS_PM_Calibration = TF_PM_HATS_L./TF_PM_HATS_R;
           
%% Finding the Calibration Transfer Functions (Mini to H.A.T.S.)
% ***NOTE*** I did not multiply the FFT by two to account for the negative
% spectrum....this is acounted for since the doubling cancels out when the
% frequency responses are divided by each other

% Averaging the respective impulse response sets

    % LEFT MINI MICROPHONE IMPULSES
    IR_MRTF_Mini_L_AVG = mean(IR_MRTF_Mini_L,2);

    % RIGHT MINI MICROPHONE IMPULSES
    IR_MRTF_Mini_R_AVG = mean(IR_MRTF_Mini_R,2);

    % RIGHT H.A.T.S. MICROPHONE IMPULSES
    IR_MRTF_HATS_L_AVG = mean(IR_MRTF_HATS_L,2);
    
    % RIGHT H.A.T.S. MICROPHONE IMPULSES
    IR_MRTF_HATS_R_AVG = mean(IR_MRTF_HATS_R,2);
    
% Finding the DFTs of the different calibration transfer functions

    % LEFT MINI MICROPHONE TRANSFER FUNCTION 
    TF_MRTF_Mini_L = fft(IR_MRTF_Mini_L_AVG, N);
    
    % RIGHT MINI MICROPHONE TRANSFER FUNCTION 
    TF_MRTF_Mini_R = fft(IR_MRTF_Mini_R_AVG, N);
    
    % LEFT H.A.T.S. MICROPHONE TRANSFER FUNCTION 
    TF_MRTF_HATS_L = fft(IR_MRTF_HATS_L_AVG, N);
    
    % RIGHT H.A.T.S. MICROPHONE TRANSFER FUNCTION 
    TF_MRTF_HATS_R = fft(IR_MRTF_HATS_R_AVG, N);
    
% Calibration transfer functions

    % LEFT SIDE (HATS_Mic_L./mini_Mic_L)
    left_MRTF_Calibration = TF_MRTF_HATS_L./TF_MRTF_Mini_L;
    
    % RIGHT SIDE (HATS_Mic_R./mini_Mic_R)
    right_MRTF_Calibration = TF_MRTF_HATS_R./TF_MRTF_Mini_R;

%% ***Temporary*** Evaluating the Quality of Average Calibration

% Finding the individual TFs for left and right microhpnes (mini & H.A.T.S.)
    
    % LEFT MINI
    mini_L_1 = fft(IR_MRTF_Mini_L(:,1),N);
    mini_L_2 = fft(IR_MRTF_Mini_L(:,2),N);
    mini_L_3 = fft(IR_MRTF_Mini_L(:,3),N);
    mini_L_4 = fft(IR_MRTF_Mini_L(:,4),N);

    % RIGHT MINI
    mini_R_1 = fft(IR_MRTF_Mini_R(:,1),N);
    mini_R_2 = fft(IR_MRTF_Mini_R(:,2),N);
    mini_R_3 = fft(IR_MRTF_Mini_R(:,3),N);
    mini_R_4 = fft(IR_MRTF_Mini_R(:,4),N);
    
    % LEFT H.A.T.S.
    HATS_L_1 = fft(IR_MRTF_HATS_L(:,1),N);
    HATS_L_2 = fft(IR_MRTF_HATS_L(:,2),N);
    HATS_L_3 = fft(IR_MRTF_HATS_L(:,3),N);
    HATS_L_4 = fft(IR_MRTF_HATS_L(:,4),N);

    % RIGHT H.A.T.S.
    HATS_R_1 = fft(IR_MRTF_HATS_R(:,1),N);
    HATS_R_2 = fft(IR_MRTF_HATS_R(:,2),N);
    HATS_R_3 = fft(IR_MRTF_HATS_R(:,3),N);
    HATS_R_4 = fft(IR_MRTF_HATS_R(:,4),N);
    
% Calibrating the individual measurements

    % LEFT CALIBRATED
    cal_mini_L_1 = mini_L_1 .* left_MRTF_Calibration;
    cal_mini_L_2 = mini_L_2 .* left_MRTF_Calibration;
    cal_mini_L_3 = mini_L_3 .* left_MRTF_Calibration;
    cal_mini_L_4 = mini_L_4 .* left_MRTF_Calibration;
    
    % RIGHT CALIBRATED
    cal_mini_R_1 = mini_R_1 .* right_MRTF_Calibration;
    cal_mini_R_2 = mini_R_2 .* right_MRTF_Calibration;
    cal_mini_R_3 = mini_R_3 .* right_MRTF_Calibration;
    cal_mini_R_4 = mini_R_4 .* right_MRTF_Calibration;

% Plotting the calibration results of the individual measurements

    % LEFT EAR
    figure
    plot(20*log10(abs(HATS_L_1(1:end/2))))
    hold on
    plot(20*log10(abs(cal_mini_L_1(1:end/2))))
    plot(20*log10(abs(HATS_L_2(1:end/2))))
    plot(20*log10(abs(cal_mini_L_2(1:end/2))))
    plot(20*log10(abs(HATS_L_3(1:end/2))))
    plot(20*log10(abs(cal_mini_L_3(1:end/2))))
    plot(20*log10(abs(HATS_L_4(1:end/2))))
    plot(20*log10(abs(cal_mini_L_4(1:end/2))))
    grid on
    title('Comparing Calibration For Left Microphones')

    % RIGHT EAR
    figure
    plot(20*log10(abs(HATS_R_1(1:end/2))))
    hold on
    plot(20*log10(abs(cal_mini_R_1(1:end/2))))
    plot(20*log10(abs(HATS_R_2(1:end/2))))
    plot(20*log10(abs(cal_mini_R_2(1:end/2))))
    plot(20*log10(abs(HATS_R_3(1:end/2))))
    plot(20*log10(abs(cal_mini_R_3(1:end/2))))
    plot(20*log10(abs(HATS_R_4(1:end/2))))
    plot(20*log10(abs(cal_mini_R_4(1:end/2))))
    grid on
    title('Comparing Calibration For Left Microphones')

%% Plotting the Calibration Transfer Functions (Phase Matching)

% Plotting phase matching calibration data

    % (mini_Mic_L./mini_Mic_R)
    figure
    subplot(2,1,1)
    plot(f,20*log10(abs(mini_PM_Calibration(1:(end/2+1)))))
    grid on
    title('(mini Mic L./mini Mic R) Magnitude (dB)')
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    subplot(2,1,2)
    plot(f,angle(mini_PM_Calibration(1:(end/2+1))))
    grid on
    title('(mini Mic L./mini Mic R) Phase (radians)')
    xlabel('Frequency (Hz)')
    ylabel('Phase (radians)')
    
    % (HATS_L./HATS_R)
    figure
    subplot(2,1,1)
    plot(f,20*log10(abs(HATS_PM_Calibration(1:(end/2+1)))))
    grid on
    title('(H.A.T.S. Mic L./H.A.T.S. Mic R) Magnitude (dB)')
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    subplot(2,1,2)
    plot(f,angle(HATS_PM_Calibration(1:(end/2+1))))
    grid on
    title('(H.A.T.S. Mic L./H.A.T.S. Mic R) Phase (radians)')
    xlabel('Frequency (Hz)')
    ylabel('Phase (radians)')
    
%% Plotting the Calibration Transfer Functions (Mini to H.A.T.S.)

    % LEFT SIDE (HATS_Mic_L./mini_Mic_L)
    figure
    subplot(2,1,1)
    plot(f,20*log10(abs(left_MRTF_Calibration(1:(end/2+1)))))
    grid on
    title('LEFT SIDE CALIBRATION (HATS Mic L./mini Mic L) Magnitude (dB)')
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    subplot(2,1,2)
    plot(f,angle(left_MRTF_Calibration(1:(end/2+1))))
    grid on
    title('LEFT SIDE CALIBRATION (HATS Mic L./mini Mic L) Phase (radians)')
    xlabel('Frequency (Hz)')
    ylabel('Phase (radians)')
    
    % RIGHT SIDE (HATS_Mic_R./mini_Mic_R)
    figure
    subplot(2,1,1)
    plot(f,20*log10(abs(right_MRTF_Calibration(1:(end/2+1)))))
    grid on
    title('RIGHT SIDE CALIBRATION (HATS Mic R./mini Mic R) Magnitude (dB)')
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    subplot(2,1,2)
    plot(f,angle(right_MRTF_Calibration(1:(end/2+1))))
    grid on
    title('RIGHT SIDE CALIBRATION (HATS Mic R./mini Mic R) Phase (radians)')
    xlabel('Frequency (Hz)')
    ylabel('Phase (radians)')
