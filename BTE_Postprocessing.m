%% BTE_Postprocessing.m
%
% This script processes the Behind the Ear (BTE) binaural room impulse 
% response data 
% 
% ** -> User Action
% *** -> Important Note (Revisit)
%
% Last updated by Michael Klunk 02/21/2021

%% Initializing the Script & (** Setting Flags **)

clear
close all
clc

% If true -> plot TFs, if flase -> no TF plots
flag_TF_Plot = false;    

% If true -> plot TFs, if flase -> no Cascade TF plots
flag_Cascade_TF = false;

%% Loading in BTE Measurement Data

n_Meas = 38;
param_Str = '/Users/michaelklunk/Programming/Audio Booth BRIR Postprocessing/Data/BTEData_2.mat';
data_Str = '/Users/michaelklunk/Documents/Msc_Dissertation_Files/Lab Day 2/';
BTE_Data = loadBRIRs(n_Meas, param_Str, data_Str);


% Removing the missed measurement, 'neg_4_pos_60_45deg.mat'
% *** This is unique to this particular set of measurements, remove for
% other scripts ***
BTE_Data(12) = [];

%% Taking the DFTs of the Measured BRIRs

% Loading in the calibration transfer functions
left_Cal = '/Users/michaelklunk/Programming/Audio Booth BRIR Postprocessing/Data/left_MRTF_Calibration.mat';
right_Cal = '/Users/michaelklunk/Programming/Audio Booth BRIR Postprocessing/Data/right_MRTF_Calibration.mat';

% Number of FFT bins
N = 2^11;

[BTE_Data,f] = processBRIRs(BTE_Data, left_Cal, right_Cal, N);

%% Deconvolution all TFs to N (Defined Above) Length Impulse Responses (See “Fast deconvolution of multichannel systems using regularization”)
% *** It may be a good idea to keep the H.A.T.S. HRTFs as full length HRTFs
% to truly encapsulate how the room would effect the CTC...Ask trever
% about this ***
% *** ALSO, I can probably stick this in a function...I don't think it
% needs to be in this script... ***

% For loop interating over all measurement sets
for index = 1:size(BTE_Data,2)
    
    % Left speaker
    BTE_Data(index).IR_Short_LEFT = ifft(BTE_Data(index).TF_LEFT, N);

    % Right speaker
    BTE_Data(index).IR_Short_RIGHT = ifft(BTE_Data(index).TF_RIGHT, N);
    
end
    
%% Using Time Least Squares Filter Design to Find the Sys Cascade Impulse Responses  (*** STARTING WITH JUST THE MINI MICS ***)

% Input Parameters

    % Sampling Frequency *** This needs to change if measurements are taken at
    % a different sample rate ***
    Fs = 48e3;

    % Size of CTC filters
    Nc = 2^11;

    % Modeling delay
    md = Nc/2;

    % Length of plant impulse responses (shortened versions)
    Nh = length(BTE_Data(1).IR_Short_LEFT(:,1));

    % Regularization filter length 
    Nr = 32;

    % Passband & stopband normalized frequencies
    f_pb = 13e3/(Fs/2);
    f_sb = 15e3/(Fs/2);

    % Maximum CTC gain (dB) Using Kabzinski Appendex: A and Masiero Appendex: A, to find the Reg Factor
    dB_PB = 30;
    dB_SB = 0;

% Preforming the cross talk cancellation and finding the sys cascade IRs
[p_11,p_12,p_21,p_22] = TD_LS_BRIR_Design(BTE_Data,Nc,Nh,md,Nr,f_pb,f_sb,dB_PB,dB_SB);

%% Finding the Cascade System Frequency Responses

% FFT Parameters
N_ = 2^16;                           % Length of FFT
dt = 1/Fs;                           % Time between samples (s)
T = dt*N_;                           % FFT Period (s)
df = Fs/N_;                          % Frequency spacing (Hz)
f_Normalized = (0:(N_/2))/N_ * 2;    % Frequencies normalized to 1 (Hz_norm)
f = f_Normalized*Fs / 2;             % Frequency (Hz)

% Finding the frequency response of all the cascade IRs
for index = 1:size(BTE_Data,2)

    P_11{index} = fft(p_11{index}, N_);
    P_12{index} = fft(p_12{index}, N_);
    P_21{index} = fft(p_21{index}, N_);
    P_22{index} = fft(p_22{index}, N_);

end

%% Evaluating CTC Preformance (Finding Applicable Metrics)
    
% Setting the lower & upper frequency limmit for preformance eval

    % Finding the index for 250 Hz
    f_low = find(f >= 250, 1);

    % Finding the index for 15,000 Hz
    f_high = find(f >= (15e3), 1);

% Length of limited frequency range
N_f_lim = length(f_low:f_high);

% Post-processing to find the preformance metrics
BTE_Data = CTC_preformanceMetrics(BTE_Data,P_11,P_12,P_21,P_22,f_low,f_high,N_f_lim,p_11,p_12,p_21,p_22,f,N_);

%% Plotting the Transfer Functions (** FLAGGED **) (** Select a Measurement **)

% ** SELECT A MEASUREMENT TO PLOT **
meas_Number = 10;

% Frequency list for the plots (Hz)
f_TF = (0:((N/2)))./N*Fs;      

if flag_TF_Plot == true
    
    figure
    semilogx(f_TF,20*log10(abs(BTE_Data(meas_Number).TF_LEFT(1:(end/2+1),1))))
    hold on
    semilogx(f_TF,20*log10(abs(BTE_Data(meas_Number).TF_LEFT(1:(end/2+1),3))))
    grid on
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    title('BTE: Left Loudspeaker to Left Ear')
    xlim([200 16000])
    legend('Mini','H.A.T.S.')

    figure
    semilogx(f_TF,20*log10(abs(BTE_Data(meas_Number).TF_LEFT(1:(end/2+1),2))))
    hold on
    semilogx(f_TF,20*log10(abs(BTE_Data(meas_Number).TF_LEFT(1:(end/2+1),3))))
    grid on
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    title('BTE: Left Loudspeaker to Right Ear')
    xlim([200 16000])
    legend('Mini','H.A.T.S.')

    figure
    semilogx(f_TF,20*log10(abs(BTE_Data(meas_Number).TF_RIGHT(1:(end/2+1),1))))
    hold on
    semilogx(f_TF,20*log10(abs(BTE_Data(meas_Number).TF_RIGHT(1:(end/2+1),3))))
    grid on
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    title('BTE: Right Loudspeaker to Left Ear')
    xlim([200 16000])
    legend('Mini','H.A.T.S.')

    figure
    semilogx(f_TF,20*log10(abs(BTE_Data(meas_Number).TF_RIGHT(1:(end/2+1),2))))
    hold on
    semilogx(f_TF,20*log10(abs(BTE_Data(meas_Number).TF_RIGHT(1:(end/2+1),4))))
    grid on
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    title('BTE: Right Loudspeaker to Right Ear')
    xlim([200 16000])
    legend('Mini','H.A.T.S.')
    
else 
    
end

%% Plotting the Cascade System Impulse & Frequency Responses **, ***
% (** FLAGGED **) 
% *** Note, you can also change from left ear to right ear ***

if flag_Cascade_TF == true
    
    figure 
    hold on
    plot(p_11{meas_Number});
    plot(p_12{meas_Number});
    plot(p_21{meas_Number});
    plot(p_22{meas_Number});
    grid on 
    title('Cascade System Impulse Responses')
    ylabel('Amplitude')
    xlabel('Sample')
    legend('p_11', 'p_12', 'p_21', 'p_22')

    figure 
    semilogx(f,20*log10(abs(P_11{meas_Number}(1:(N_/2+1)))));
    hold on 
    semilogx(f,20*log10(abs(P_21{meas_Number}(1:(N_/2+1)))));
    xlim([100 15000])
    grid on

end
