%% BTE_Postprocessing.m
%
% This script processing the Behind the Ear (BTE) impulse response data 
%
% Last updated by Michael Klunk 02/16/2021

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
N = 2^8;                    % Number of FFT Points
dt = 1/Fs;                  % Delta t (s)
T = dt*N;                   % Sampling Period (s)
df = 1/T;                   % Delta f (Hz)
f = (0:((N/2)-0))./N*Fs;    % Frequency vector (Hz)

% Loading in the calibration transfer functions
load('left_MRTF_Calibration.mat');
load('right_MRTF_Calibration.mat');

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
    concat_TF_BTE_Right_Speaker_TEST = [right_Speaker_Left_Ear, right_Speaker_Right_Ear, TF_BTE_Right_SPEAKER_TEST(:,3:4)];
    
%% Plotting the Transfer Functions

figure
plot(f,20*log10(abs(concat_TF_BTE_Left_SPEAKER_TEST((1:(end/2+1)),1))))
hold on
plot(f,20*log10(abs(concat_TF_BTE_Left_SPEAKER_TEST((1:(end/2+1)),3))))
grid on
legend('Mini','H.A.T.S.')

figure
plot(f,20*log10(abs(concat_TF_BTE_Left_SPEAKER_TEST((1:(end/2+1)),2))))
hold on
plot(f,20*log10(abs(concat_TF_BTE_Left_SPEAKER_TEST((1:(end/2+1)),4))))
grid on
legend('Mini','H.A.T.S.')

figure
plot(f,20*log10(abs(concat_TF_BTE_Right_Speaker_TEST((1:(end/2+1)),1))))
hold on
plot(f,20*log10(abs(concat_TF_BTE_Right_Speaker_TEST((1:(end/2+1)),3))))
grid on
legend('Mini','H.A.T.S.')

figure
plot(f,20*log10(abs(concat_TF_BTE_Right_Speaker_TEST((1:(end/2+1)),2))))
hold on
plot(f,20*log10(abs(concat_TF_BTE_Right_Speaker_TEST((1:(end/2+1)),4))))
grid on
legend('Mini','H.A.T.S.')