function [BRIR_Data] = CTC_preformanceMetrics(BRIR_Data,P_11,P_12,P_21,P_22,f_low,f_high,N_f_lim,p_11,p_12,p_21,p_22,f,N_)
% CTC_preformanceMetrics
%
%   Input parameters:
% 
%       P_11,P_12,P_21,P_22:     Cascade System Frequency Responses
%
%       f_low,f_high:            Lower & upper frequency limmit for preformance eval
%      
%       BRIR_Data:               Struct containing the following fields:
%                                   - abscissa 
%                                   - ordinate
%                                   - rotation
%                                   - IR_LEFT (IRs from the left speaker)
%                                   - IR_RIGHT (IRs from the right speaker)
%                                   - TF_LEFT (TFs from the left speaker)
%                                   - TF_RIGHT (TFs from the right speaker)
%                                   - IR_Short_LEFT (Shortened IRs from the
%                                       left speaker)                           
%                                   - IR_Short_RIGHT (Shortened IRs from the
%                                       right speaker)
%
%       N_f_lim:                 Number of Bins in evaluation spectrum
%
%       p_11,p_12,p_21,p_22:     Cascade System Impulse Responses
%
%       f:                       Frequency List (Hz)
%
%       N_:                      Number of FFT Points
%
%   Output parameters:
% 
%       BRIR_Data:               Struct containing the following fields:                             
% 
% #Author: Michael Klunk 
% #Date: Monday, February 21st, 2022

%% Maximum Deviation

for index = 1:size(BRIR_Data,2)
    
    % Calculating the max deviation for desired direct path cascade systems
    max_Dev_11 = 20*log10(max(abs(P_11{index}(f_low:f_high,1))));
    max_Dev_22 = 20*log10(max(abs(P_22{index}(f_low:f_high,1))));
    BRIR_Data(index).max_Dev = max([max_Dev_11, max_Dev_22]);
    
end

%% Maximum Crosstalk Contribution 

for index = 1:size(BRIR_Data,2)
    
    % Calculating the max cross-talk contribution for the applicable cascade systems
    max_CT_12 = 20*log10(max(abs(P_12{index}(f_low:f_high,1))));
    max_CT_21 = 20*log10(max(abs(P_21{index}(f_low:f_high,1))));
    BRIR_Data(index).max_CT = max([max_CT_12, max_CT_21]);
    
end

%% Minimum Channel Seperation

for index = 1:size(BRIR_Data,2)
    
    % Calculating the min channel seperation for the applicable cascade systems
    min_CS_11 = 20*log10(min(abs(P_11{index}(f_low:f_high,1) ./ P_12{index}(f_low:f_high,1))));                      
    min_CS_22 = 20*log10(min(abs(P_22{index}(f_low:f_high,1) ./ P_21{index}(f_low:f_high,1))));
    BRIR_Data(index).min_CS = min([min_CS_11, min_CS_22]);
    
end

%% Mean Channel Seperation

for index = 1:size(BRIR_Data,2)

    % Calculating the mean channel seperation
    BRIR_Data(index).mean_CS = 20/(N_f_lim*2*(2-1)) * (sum(log10(abs(P_11{index}(f_low:f_high,1)...
        ./ P_12{index}(f_low:f_high,1)))) + sum(log10(abs(P_22{index}(f_low:f_high,1)...
        ./ P_21{index}(f_low:f_high,1)))));

end

%% Maximum Interaural Phase Difference Distortion *** This is a bit messy, Clean it Up ***

for index = 1:size(BRIR_Data,2)
    
    P_1112 = fft((p_11{index}+p_12{index}), N_); P_1112 = P_1112(1:(end/2+1)); 
    phase_Delay_P_1112 = -unwrap((angle(P_1112)))./(f'*2*pi);

    P_2122 = fft((p_21{index}+p_22{index}), N_); P_2122 = P_2122(1:(end/2+1));
    phase_Delay_P_2122 = -unwrap((angle(P_2122)))./(f'*2*pi);

    BRIR_Data(index).max_IPDD = max(abs(phase_Delay_P_1112(f_low:f_high) - phase_Delay_P_2122(f_low:f_high)));
    
end




