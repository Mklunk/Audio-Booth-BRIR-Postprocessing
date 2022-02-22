function [p_11,p_12,p_21,p_22] = TD_LS_BRIR_Design(BRIR_Data,Nc,Nh,md,Nr,f_pb,f_sb,dB_PB,dB_SB)
% TD_LS_BRIR_Design
%
%   Input parameters:
%
%       BRIR_Data:          Struct containing the following fields:
%                               - abscissa 
%                               - ordinate
%                               - rotation
%                               - IR_LEFT (IRs from the left speaker)
%                               - IR_RIGHT (IRs from the right speaker)
%                               - TF_LEFT (TFs from the left speaker)
%                               - TF_RIGHT (TFs from the right speaker)
%                               - IR_Short_LEFT (Shortened IRs from the
%                                   left speaker)                           
%                               - IR_Short_RIGHT (Shortened IRs from the
%                                   right speaker)
%          
%       Nc:                 Size of CTC filters
%
%       Nh:                 Length of plant impulse responses (shortened
%                           versions)
%
%       md:                 Modeling delay to ensure causality
%
%       Nr:                 Regularization filter length
%
%       f_pb:               Passband normalized frequency
%
%       f_sb:               Stopband normalized frequency
%
%       dB_PB:              Maximum passband gain (dB)
%
%       dB_SB:              Maximum stopband gain (dB)
% 
%   Output parameters:
% 
%       p_11,p_12,p_21,p_22: System Cascade Impulse Responses 
%                            (where p_12 -> receiver 1, channel 2)
% 
% #Author: Michael Klunk 
% #Date: Monday, February 21st, 2022

% Creating the convolution matrices for the BRTF systems 
% (***STARTING WITH JUST THE MINI MICS***)
% (LE_LL) -> (Left Ear, Left Speaker)
for index = 1:size(BRIR_Data,2)
    
    h_LE_LS_Conv{index} = convmtx(BRIR_Data(index).IR_Short_LEFT(:,1), Nc);
    h_RE_LS_Conv{index} = convmtx(BRIR_Data(index).IR_Short_LEFT(:,2), Nc);
    h_LE_RS_Conv{index} = convmtx(BRIR_Data(index).IR_Short_RIGHT(:,1), Nc);
    h_RE_RS_Conv{index} = convmtx(BRIR_Data(index).IR_Short_RIGHT(:,2), Nc);
    
    % Composite convolution matrix 
    h_conv_Comp{index} = [h_LE_LS_Conv{index}, h_LE_RS_Conv{index};...
        h_RE_LS_Conv{index}, h_RE_RS_Conv{index}];

end

%% Creating the Desired System Response Vectors 

% Finding length of the desired response
Nd = Nh + Nc - 1;

% Creating the desired responses, d (where d_rec1_chan2 -> receiver 1, channel 2)
d_rec1_chan1 = zeros(Nd,1); d_rec1_chan1(md+1,1) = 1;
d_rec1_chan2 = zeros(Nd,1);
d_rec2_chan1 = zeros(Nd,1);
d_rec2_chan2 = zeros(Nd,1); d_rec2_chan2(md+1,1) = 1;

% Creating the cascaded desired responses for input 1 & 2
d_input1 = [d_rec1_chan1; d_rec2_chan1];
d_input2 = [d_rec1_chan2; d_rec2_chan2];

%% Creating the Regularization Filter 

% Using Kabzinski Appendex: A and Masiero Appendex: A, to find the Reg Factor
a_PB = 1/((2*10^(dB_PB/20))^2);
a_SB = 1/((2*10^(dB_SB/20))^2);

% Creating the filter using firls
b = firls(Nr, [0 f_pb f_sb 1], [a_PB a_PB a_SB a_SB]); b = b.';

% Creating the convolution matrix for filter, b
b_conv = convmtx(b, Nc);

% Leveraging the kronecker product
b_Kron = kron(eye(2),(transpose(b_conv)*b_conv));

%% Solving for the Linear Least-Squares Filter
% ** Renaming the filters I think might be a good idea ***

for index = 1:size(BRIR_Data,2)
    
    % SOLVING FOR INPUT ONE
        % A from -> A*x = B 
        A_linEq_1 = transpose(h_conv_Comp{index})*h_conv_Comp{index} + b_Kron;

        % B from -> A*x = B 
        B_linEq_1 = transpose(h_conv_Comp{index})*d_input1;

        % Solving the sys of lin eqs
        c_1{index} = A_linEq_1\B_linEq_1;

        % Spliting c_1 into c_11 & c_21
        c_11{index} = c_1{index}(1:floor(end/2));
        c_21{index} = c_1{index}(floor(end/2)+1:end);

    % SOLVING FOR INPUT TWO
        % A from -> A*x = B 
        A_linEq_2 = transpose(h_conv_Comp{index})*h_conv_Comp{index} + b_Kron;

        % B from -> A*x = B 
        B_linEq_2 = transpose(h_conv_Comp{index})*d_input2;

        % Solving the sys of lin eqs
        c_2{index} = A_linEq_2\B_linEq_2;

        % Spliting c_2 into c_12 & c_22
        c_12{index} = c_2{index}(1:floor(end/2));
        c_22{index} = c_2{index}(floor(end/2)+1:end);

end

%% Finding the Cascade Impulse Responses

for index = 1:size(BRIR_Data,2)
    
    % Creating the cascade responses, p (where p_12 -> receiver 1, channel 2)
    p_11{index} = conv(BRIR_Data(index).IR_Short_LEFT(:,3), c_11{index}) + ...
        conv(BRIR_Data(index).IR_Short_RIGHT(:,3), c_21{index});
    p_12{index} = conv(BRIR_Data(index).IR_Short_LEFT(:,3), c_12{index}) + ...
        conv(BRIR_Data(index).IR_Short_RIGHT(:,3), c_22{index});
    p_21{index} = conv(BRIR_Data(index).IR_Short_LEFT(:,4), c_11{index}) + ...
        conv(BRIR_Data(index).IR_Short_RIGHT(:,4), c_21{index});
    p_22{index} = conv(BRIR_Data(index).IR_Short_LEFT(:,4), c_12{index}) + ...
        conv(BRIR_Data(index).IR_Short_RIGHT(:,4), c_22{index});
    
end


