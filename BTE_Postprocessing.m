%% BTE_Postprocessing.m
%
% This script processes the Behind the Ear (BTE) impulse response data 
%
% Last updated by Michael Klunk 02/17/2021

%% Initializing the Script

clear
close all
clc

%% Loading in BTE Measurement Data

% Number of Measurements Taken
num_Measurements = 40;

% Creating empty cell array to initialize BTE struct
struct_Init = cell(1,num_Measurements);
for index=1:num_Measurements
   struct_Init{index} = 0;
end

% Inititalizing fields and values for BTE data
field1 = 'abscissa';  value1 = struct_Init;
field2 = 'ordinate';  value2 = struct_Init;
field3 = 'rotation';  value3 = struct_Init;
field4 = 'IR_LEFT';  value4 = {zeros(48000,4)};
field5 = 'IR_RIGHT';  value5 = {zeros(48000,4)};
field6 = 'measurement_Iteration'; value6 = 1;

% Inititalizing BTE struct
BTE_Data = struct(field1,value1,field2,value2,field3,value3,field4,...
    value4,field5,value5,field6,value6);

% Clearing unneeded initiazation variables
clear field1 field2 field3 field4 field5 field6 value1 value2 value3...
    value4 value5 value6 struct_Init

% Loading in raw BTE Data as a cell array
load('BTEData.mat');

% For loop loading all of the raw BTE data into the struct
for index = 1:38
    
    try
        data_TEMP = load(strcat('/Users/michaelklunk/Documents/Msc_Dissertation_Files/Lab Day 2/',BTEData1{index,4},'.mat'));
    catch
        disp(strcat("The file ",BTEData1{index,4},'.mat',' does not exist'));
        continue
    end
    
    BTE_Data(index).abscissa = BTEData1{index,1};
    BTE_Data(index).ordinate = BTEData1{index,2};
    BTE_Data(index).rotation = BTEData1{index,3};
    BTE_Data(index).IR_LEFT = data_TEMP.data(1).IR;
    BTE_Data(index).IR_RIGHT = data_TEMP.data(2).IR;
end

% Removing the missed measurement, 'neg_4_pos_60_45deg.mat'
BTE_Data(12) = [];

% Manually adding the few second measurement iterations

    % LOADING 'neg_4_0_0deg_2.mat'
    data_TEMP = load('/Users/michaelklunk/Documents/Msc_Dissertation_Files/Lab Day 2/neg_4_0_0deg_2.mat');
    BTE_Data(38).abscissa = 1.4635;
    BTE_Data(38).ordinate = 1.4625;
    BTE_Data(38).rotation = 0;
    BTE_Data(38).IR_LEFT = data_TEMP.data(1).IR;
    BTE_Data(38).IR_RIGHT = data_TEMP.data(2).IR;
    BTE_Data(38).measurement_Iteration = 2;
    
    % LOADING 'neg_4_pos_60_90deg_2.mat'
    data_TEMP = load('/Users/michaelklunk/Documents/Msc_Dissertation_Files/Lab Day 2/neg_4_pos_60_90deg_2.mat');
    BTE_Data(39).abscissa = 1.4635;
    BTE_Data(39).ordinate = 2.0625;
    BTE_Data(39).rotation = 90;
    BTE_Data(39).IR_LEFT = data_TEMP.data(1).IR;
    BTE_Data(39).IR_RIGHT = data_TEMP.data(2).IR;
    BTE_Data(39).measurement_Iteration = 2;
    
    % LOADING 'neg_4_pos_60_270deg_2.mat'
    data_TEMP = load('/Users/michaelklunk/Documents/Msc_Dissertation_Files/Lab Day 2/neg_4_pos_60_270deg_2.mat');
    BTE_Data(40).abscissa = 1.4635;
    BTE_Data(40).ordinate = 2.0625;
    BTE_Data(40).rotation = 270;
    BTE_Data(40).IR_LEFT = data_TEMP.data(1).IR;
    BTE_Data(40).IR_RIGHT = data_TEMP.data(2).IR;
    BTE_Data(40).measurement_Iteration = 2;
    
% Clearing unneed imported data variables
clear BTEData1 data_TEMP

%% Taking the DFT of a Measurement (***THIS IS A TEST***)
% ***NOTE*** Only using N = 2^8 for now...if this changes I will need to
% also bring in higher resolution calibration transfer functions

% FFT Parameters
Fs = 48000;                 % Sampling frequency (Hz)
N = 2^11;                   % Number of FFT Points
dt = 1/Fs;                  % Delta t (s)
T = dt*N;                   % Sampling Period (s)
df = 1/T;                   % Delta f (Hz)
f = (0:((N/2)-0))./N*Fs;    % Frequency vector (Hz)
t = dt*(0:(N-1));           % Time Vector (s)

% Loading in the calibration transfer functions
load('left_MRTF_Calibration.mat');
load('right_MRTF_Calibration.mat');

% Converting the first measurement into the frequency domain
TF_BTE_Left_SPEAKER_TEST = fft(BTE_Data(40).IR_LEFT,N);
TF_BTE_Right_SPEAKER_TEST = fft(BTE_Data(40).IR_RIGHT,N);

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

%% Deconvolution all TFs to N (2^11) Length Impulse Responses (See “Fast deconvolution ofmultichannel systems using regularization”)
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
Nc = 2048;

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
md = 256;

% Creating the desired responses, d (where d_rec1_chan2 -> receiver 1, channel 2)
d_rec1_chan1 = zeros(Nd,1); d_rec1_chan1(md+1,1) = 1;
d_rec1_chan2 = zeros(Nd,1);
d_rec2_chan1 = zeros(Nd,1);
d_22 = zeros(Nd,1); d_22(md+1,1) = 1;

% Creating the cascaded desired responses for input 1 & 2
d_input1 = [d_rec1_chan1; d_rec2_chan1];
d_input2 = [d_rec1_chan2; d_22];

%% Creating the Regularization Filter 

% Regularization filter length 
Nr = 32;

% Passband & stopband normalized frequencies
f_pb = 13e3/(Fs/2);
f_sb = 13.25e3/(Fs/2);

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
p_11 = conv(concat_IR_BTE_Left_SPEAKER_TEST(:,3), c_11) + ...
    conv(concat_IR_BTE_Right_SPEAKER_TEST(:,3), c_21);
p_12 = conv(concat_IR_BTE_Left_SPEAKER_TEST(:,3), c_12) + ...
    conv(concat_IR_BTE_Right_SPEAKER_TEST(:,3), c_22);
p_21 = conv(concat_IR_BTE_Left_SPEAKER_TEST(:,4), c_11) + ...
    conv(concat_IR_BTE_Right_SPEAKER_TEST(:,4), c_21);
p_22 = conv(concat_IR_BTE_Left_SPEAKER_TEST(:,4), c_12) + ...
    conv(concat_IR_BTE_Right_SPEAKER_TEST(:,4), c_22);

%% Calculating Cascade System Frequency Responses

% FFT Parameters
N = 2^12;                           % Length of FFT
dt = 1/Fs;                          % Time between samples (s)
T = dt*N;                           % FFT Period (s)
df = Fs/N;                          % Frequency spacing (Hz)
f_Normalized = (0:(N/2+1))/N * 2;   % Frequencies normalized to 1 (Hz_norm)
f = f_Normalized*Fs / 2;            % Frequency (Hz)

% Finding the frequency response of all the cascade IRs
P_11 = fft(p_11, N);
P_12 = fft(p_12, N);
P_21 = fft(p_21, N);
P_22 = fft(p_22, N);

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
hold on 
plot(f,20*log10(abs(P_11(1:(N/2+2)))));
plot(f,20*log10(abs(P_21(1:(N/2+2)))));
xlim([200 16000])
grid on


