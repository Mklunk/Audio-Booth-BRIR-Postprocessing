%% BTE_Postprocessing.m
%
% This script processes the Behind the Ear (BTE) binaural room impulse 
% response data 
%
% Last updated by Michael Klunk 02/21/2021

%% Initializing the Script

clear
close all
clc

%% Loading in BTE Measurement Data

n_Meas = 38;
param_Str = '/Users/michaelklunk/Programming/Audio Booth BRIR Postprocessing/Data/BTEData_2.mat';
data_Str = '/Users/michaelklunk/Documents/Msc_Dissertation_Files/Lab Day 2/';
BTE_Data = loadBRIRs(n_Meas, param_Str, data_Str);

% Removing the missed measurement, 'neg_4_pos_60_45deg.mat'
BTE_Data(12) = [];

%% Taking the DFT of a Measurement (***THIS IS A TEST***)
% ***NOTE*** Only using N = 2^8 for now...if this changes I will need to
% also bring in higher resolution calibration transfer functions

% FFT Parameters
Fs = 48000;                 % Sampling frequency (Hz)
N = 2^11;                   % Number of FFT Points
dt = 1/Fs;                  % Delta t (s)
T = dt*N;                   % Sampling Period (s)
df = 1/T;                   % Delta f (Hz)
f = (0:((N/2)))./N*Fs;    % Frequency vector (Hz)
t = dt*(0:(N-1));           % Time Vector (s)

% Loading in the calibration transfer functions
load('/Users/michaelklunk/Programming/Audio Booth BRIR Postprocessing/Data/left_MRTF_Calibration.mat');
load('/Users/michaelklunk/Programming/Audio Booth BRIR Postprocessing/Data/right_MRTF_Calibration.mat');

% WINDOWING OPTION
%                 % Adding a half hamming window to fade out the last 12 % of the samples to
%                 % avoid leakage? 
%                     % Creating the Window (Make a multiple of 0.2) <- you might not have to!
%                     w_PERCENT = 0.95;
%                     w = hann(ceil(N*w_PERCENT*2));
%                     w = [ones(N-ceil(N*w_PERCENT),1);w(end/2+1:end)];
%                 % Converting the first measurement into the frequency domain
%                     TF_BTE_Left_SPEAKER_TEST = fft(BTE_Data(1).IR_LEFT(1:N,:).*w,N);
%                     TF_BTE_Right_SPEAKER_TEST = fft(BTE_Data(1).IR_RIGHT(1:N,:).*w,N);

% NON-WINDOWING OPTION
% Converting the first measurement into the frequency domain
    TF_BTE_Left_SPEAKER_TEST = fft(BTE_Data(1).IR_LEFT,N);
    TF_BTE_Right_SPEAKER_TEST = fft(BTE_Data(1).IR_RIGHT,N);
                
% Calibrating the mini microphone data

    % LEFT EAR
     left_Speaker_Left_Ear = TF_BTE_Left_SPEAKER_TEST(:,1).*left_MRTF_Calibration;
     right_Speaker_Left_Ear = TF_BTE_Right_SPEAKER_TEST(:,1).*left_MRTF_Calibration;
     
    % RIGHT EAR
     left_Speaker_Right_Ear= TF_BTE_Left_SPEAKER_TEST(:,2).*right_MRTF_Calibration;
     right_Speaker_Right_Ear = TF_BTE_Right_SPEAKER_TEST(:,2).*right_MRTF_Calibration;
     
% Concatenating calibrated mini microphone data with HATS data

    % LEFT SPEAKER
    concat_TF_BTE_Left_SPEAKER_TEST = [left_Speaker_Left_Ear, left_Speaker_Right_Ear, TF_BTE_Left_SPEAKER_TEST(:,3:4)];
    concat_TF_BTE_Right_SPEAKER_TEST = [right_Speaker_Left_Ear, right_Speaker_Right_Ear, TF_BTE_Right_SPEAKER_TEST(:,3:4)];
    
%% Plotting the Transfer Functions

figure
semilogx(f,20*log10(abs(concat_TF_BTE_Left_SPEAKER_TEST((1:(end/2+1)),1))))
hold on
semilogx(f,20*log10(abs(concat_TF_BTE_Left_SPEAKER_TEST((1:(end/2+1)),3))))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('BTE: Left Loudspeaker to Left Ear')
xlim([200 16000])
legend('Mini','H.A.T.S.')

figure
semilogx(f,20*log10(abs(concat_TF_BTE_Left_SPEAKER_TEST((1:(end/2+1)),2))))
hold on
semilogx(f,20*log10(abs(concat_TF_BTE_Left_SPEAKER_TEST((1:(end/2+1)),4))))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('BTE: Left Loudspeaker to Right Ear')
xlim([200 16000])
legend('Mini','H.A.T.S.')

figure
semilogx(f,20*log10(abs(concat_TF_BTE_Right_SPEAKER_TEST((1:(end/2+1)),1))))
hold on
semilogx(f,20*log10(abs(concat_TF_BTE_Right_SPEAKER_TEST((1:(end/2+1)),3))))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('BTE: Right Loudspeaker to Left Ear')
xlim([200 16000])
legend('Mini','H.A.T.S.')

figure
semilogx(f,20*log10(abs(concat_TF_BTE_Right_SPEAKER_TEST((1:(end/2+1)),2))))
hold on
semilogx(f,20*log10(abs(concat_TF_BTE_Right_SPEAKER_TEST((1:(end/2+1)),4))))
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('BTE: Right Loudspeaker to Right Ear')
xlim([200 16000])
legend('Mini','H.A.T.S.')

%% Deconvolution all TFs to N (2^11) Length Impulse Responses (See “Fast deconvolution of multichannel systems using regularization”)
% ***It may be a good idea to keep the H.A.T.S. HRTFs as full length HRTFs
% to truly encapsulate how the room would effect the CTC...Ask trever
% about this

% Left speaker
concat_IR_BTE_Left_SPEAKER_TEST = ifft(concat_TF_BTE_Left_SPEAKER_TEST, N);

% Right speaker
concat_IR_BTE_Right_SPEAKER_TEST = ifft(concat_TF_BTE_Right_SPEAKER_TEST, N);

%% Creating the Convolution Matrices for the HRTF Systems (***STARTING WITH JUST THE MINI MICS***)
% (LE_LL) -> (Left Ear, Left Speaker)

% Size of CTC filters
Nc = 2^11;

% Size of Impulse Responses 
Nh = length(concat_IR_BTE_Left_SPEAKER_TEST(:,1));

h_LE_LS_Conv = convmtx(concat_IR_BTE_Left_SPEAKER_TEST(:,1), Nc);
h_RE_LS_Conv = convmtx(concat_IR_BTE_Left_SPEAKER_TEST(:,2), Nc);
h_LE_RS_Conv = convmtx(concat_IR_BTE_Right_SPEAKER_TEST(:,1), Nc);
h_RE_RS_Conv = convmtx(concat_IR_BTE_Right_SPEAKER_TEST(:,2), Nc);

% Composite convolution matrix
h_conv_Comp = [h_LE_LS_Conv, h_LE_RS_Conv; h_RE_LS_Conv, h_RE_RS_Conv];

%% Creating the Desired System Response Vectors 

% Finding length of the desired response
Nd = Nh + Nc - 1;

% Modeling delay
md = Nc/2;

% Creating the desired responses, d (where d_rec1_chan2 -> receiver 1, channel 2)
d_rec1_chan1 = zeros(Nd,1); d_rec1_chan1(md+1,1) = 1;
d_rec1_chan2 = zeros(Nd,1);
d_rec2_chan1 = zeros(Nd,1);
d_rec2_chan2 = zeros(Nd,1); d_rec2_chan2(md+1,1) = 1;

% Creating the cascaded desired responses for input 1 & 2
d_input1 = [d_rec1_chan1; d_rec2_chan1];
d_input2 = [d_rec1_chan2; d_rec2_chan2];

%% Creating the Regularization Filter 

% Regularization filter length 
Nr = 32;

% Passband & stopband normalized frequencies
f_pb = 13e3/(Fs/2);
f_sb = 15e3/(Fs/2);

% Using Kabzinski Appendex: A and Masiero Appendex: A, to find the Reg Factor
dB_PB = 30;
dB_SB = 0;
a_PB = 1/((2*10^(dB_PB/20))^2);
a_SB = 1/((2*10^(dB_SB/20))^2);

% Creating the filter using firls
b = firls(Nr, [0 f_pb f_sb 1], [a_PB a_PB a_SB a_SB]); b = b.';

% Creating the convolution matrix for filter, b
b_conv = convmtx(b, Nc);

% Leveraging the kronecker product
b_Kron = kron(eye(2),(transpose(b_conv)*b_conv));

%% Solving for the Linear Least-Squares Filter
% **Renaming the filters I think might be a good idea***

% SOLVING FOR INPUT ONE
    % A from -> A*x = B 
    A_linEq_1 = transpose(h_conv_Comp)*h_conv_Comp + b_Kron;

    % B from -> A*x = B 
    B_linEq_1 = transpose(h_conv_Comp)*d_input1;

    % Solving the sys of lin eqs
    c_1 = A_linEq_1\B_linEq_1;
    
    % Spliting c_1 into c_11 & c_21
    c_11 = c_1(1:floor(end/2));
    c_21 = c_1(floor(end/2)+1:end);
    
% SOLVING FOR INPUT TWO
    % A from -> A*x = B 
    A_linEq_2 = transpose(h_conv_Comp)*h_conv_Comp + b_Kron;

    % B from -> A*x = B 
    B_linEq_2 = transpose(h_conv_Comp)*d_input2;

    % Solving the sys of lin eqs
    c_2 = A_linEq_2\B_linEq_2;
    
    % Spliting c_2 into c_12 & c_22
    c_12 = c_2(1:floor(end/2));
    c_22 = c_2(floor(end/2)+1:end);
    
%% Finding the Cascade Impulse Responses

% Creating the cascade responses, p (where p_12 -> receiver 1, channel 2)
p_11 = conv(concat_IR_BTE_Left_SPEAKER_TEST(:,1), c_11) + ...
    conv(concat_IR_BTE_Right_SPEAKER_TEST(:,1), c_21);
p_12 = conv(concat_IR_BTE_Left_SPEAKER_TEST(:,1), c_12) + ...
    conv(concat_IR_BTE_Right_SPEAKER_TEST(:,1), c_22);
p_21 = conv(concat_IR_BTE_Left_SPEAKER_TEST(:,2), c_11) + ...
    conv(concat_IR_BTE_Right_SPEAKER_TEST(:,2), c_21);
p_22 = conv(concat_IR_BTE_Left_SPEAKER_TEST(:,2), c_12) + ...
    conv(concat_IR_BTE_Right_SPEAKER_TEST(:,2), c_22);

%% Calculating Cascade System Frequency Responses

% FFT Parameters
N_ = 2^16;                           % Length of FFT
dt = 1/Fs;                          % Time between samples (s)
T = dt*N_;                           % FFT Period (s)
df = Fs/N_;                          % Frequency spacing (Hz)
f_Normalized = (0:(N_/2))/N_ * 2;   % Frequencies normalized to 1 (Hz_norm)
f = f_Normalized*Fs / 2;            % Frequency (Hz)

% Finding the frequency response of all the cascade IRs
P_11 = fft(p_11, N_);
P_12 = fft(p_12, N_);
P_21 = fft(p_21, N_);
P_22 = fft(p_22, N_);

%% Plotting the Cascade System Impulse & Frequency Responses

figure 
hold on
plot(p_11);
plot(p_12);
plot(p_21);
plot(p_22);
grid on 
title('Cascade System Impulse Responses')
ylabel('Amplitude')
xlabel('Sample')
legend('p_11', 'p_12', 'p_21', 'p_22')

figure 
semilogx(f,20*log10(abs(P_11(1:(N_/2+1)))));
hold on 
semilogx(f,20*log10(abs(P_21(1:(N_/2+1)))));
xlim([100 15000])
grid on

%% Maximum Deviation 
    
% Finding the index for 250 Hz
f_250 = find(f >= 250, 1);

% Finding the index for 15,000 Hz
f_15k = find(f >= (5e3), 1);

% Length of limited frequency range
N_f_lim = length(f_250:f_15k);

% Calculating the max deviation for the applicable cascade systems
max_Dev_11 = 20*log10(max(abs(P_11(f_250:f_15k,1))));
max_Dev_22 = 20*log10(max(abs(P_22(f_250:f_15k,1))));
max_Dev = max([max_Dev_11, max_Dev_22]);
    
%% Maximum Crosstalk Contribution 

% Calculating the max cross-talk contribution for the applicable cascade systems
max_CT_12 = 20*log10(max(abs(P_12(f_250:f_15k,1))));
max_CT_21 = 20*log10(max(abs(P_21(f_250:f_15k,1))));
max_CT = max([max_CT_12, max_CT_21]);

%% Minimum Channel Seperation

% Calculating the min channel seperation for the applicable cascade systems
min_CS_11 = 20*log10(min(abs(P_11(f_250:f_15k,1) ./ P_12(f_250:f_15k,1))));                      
min_CS_22 = 20*log10(min(abs(P_22(f_250:f_15k,1) ./ P_21(f_250:f_15k,1))));

%% Mean Channel Seperation

mean_CS = 20/(N_f_lim*2*(2-1)) * (sum(log10(abs(P_11(f_250:f_15k,1)...
    ./ P_12(f_250:f_15k,1)))) + sum(log10(abs(P_22(f_250:f_15k,1)...
    ./ P_21(f_250:f_15k,1)))));

%% Maximum Interaural Phase Difference Distortion

P_1112 = fft((p_11+p_12), N_); P_1112 = P_1112(1:(end/2+1)); 
phase_Delay_P_1112 = -angle(P_1112)./(f'*2*pi);


P_2122 = fft((p_21+p_22), N_); P_2122 = P_2122(1:(end/2+1));
phase_Delay_P_2122 = -angle(P_2122)./(f'*2*pi);


max_IPDD = max(abs(phase_Delay_P_1112(700:end) - phase_Delay_P_2122(700:end)));


