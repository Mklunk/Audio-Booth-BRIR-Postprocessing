function [BRIR_Data, f] = processBRIRs(BRIR_Data, left_Cal, right_Cal, N)
% processBRIRs
%
% This function processes the BRIRs passed into the function (including calibraiton).  
%
%   Input parameters:
%
%       BRIR_Data:     Struct containing the following fields:
%                               - abscissa 
%                               - ordinate
%                               - rotation
%                               - IR_LEFT (IRs from the left speaker)
%                               - IR_RIGHT (IRs from the right speaker)
% 
%       left_Cal:           Calibration Transfer function for left mini
%                           microphone
%
%       right_Cal:          File name of the different measurement
%                           positions and file names
% 
%   Output parameters:
%
%       BRIR_Data:     Struct containing the following fields:
%                               - abscissa 
%                               - ordinate
%                               - rotation
%                               - IR_LEFT (IRs from the left speaker)
%                               - IR_RIGHT (IRs from the right speaker)
%                               - TF_LEFT (TFs from the left speaker)
%                               - TF_RIGHT (TFs from the right speaker)
% 
% #Author: Michael Klunk 
% #Date: Monday, February 21st, 2022

% FFT Parameters
Fs = 48000;                 % Sampling frequency (Hz)
dt = 1/Fs;                  % Delta t (s)
T = dt*N;                   % Sampling Period (s)
df = 1/T;                   % Delta f (Hz)
f = (0:((N/2)))./N*Fs;      % Frequency vector (Hz)
t = dt*(0:(N-1));           % Time Vector (s)

% Loading in the calibration transfer functions
load(left_Cal);
load(right_Cal);

% ***POTENTIALLY ADD WINDOWING FUNCTION HERE, It is commented out below***


% Taking the FFTs of the BRIRs and Loading them into 
for index = 1:size(BRIR_Data,2)
    
    BRIR_Data(index).TF_LEFT = fft(BRIR_Data(index).IR_LEFT,N);
    BRIR_Data(index).TF_RIGHT = fft(BRIR_Data(index).IR_RIGHT,N);
    
    % Calibrating the mini microphone data

        % LEFT EAR
        BRIR_Data(index).TF_LEFT(:,1) = BRIR_Data(index).TF_LEFT(:,1).*left_MRTF_Calibration;
        BRIR_Data(index).TF_RIGHT(:,1) = BRIR_Data(index).TF_RIGHT(:,1).*left_MRTF_Calibration;
     
        % RIGHT EAR
        BRIR_Data(index).TF_LEFT(:,2) = BRIR_Data(index).TF_LEFT(:,2).*right_MRTF_Calibration;
        BRIR_Data(index).TF_RIGHT(:,2) = BRIR_Data(index).TF_RIGHT(:,2).*right_MRTF_Calibration;

end







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
