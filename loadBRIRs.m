function [BRIR_Data] = loadBRIRs(n_Meas, param_Str, data_Str)
% loadBRIRs
%
%   Input parameters:
% 
%       n_Measurements:     Number of measurements taken (four BRIRs per set)
%
%       param_Str:          File name of the different measurement
%                           positions and file names
%
%       data_Str:           File name of the .mat file containing the
%                           binaural room impulse response (BRIR) data
% 
%   Output parameters:
% 
%       BRIR_Data:          Struct containing the following fields:
%                               - abscissa 
%                               - ordinate
%                               - rotation
%                               - IR_LEFT (IRs from the left speaker)
%                               - IR_RIGHT (IRs from the right speaker)
% 
% #Author: Michael Klunk 
% #Date: Monday, February 21st, 2022

% Creating empty cell array to initialize BRIR struct
struct_Init = cell(1,n_Meas);
for index=1:n_Meas
   struct_Init{index} = 0;
end

% Inititalizing fields and values for BRIR data
field1 = 'abscissa';  value1 = struct_Init;
field2 = 'ordinate';  value2 = struct_Init;
field3 = 'rotation';  value3 = struct_Init;
field4 = 'IR_LEFT';  value4 = {zeros(48000,4)};
field5 = 'IR_RIGHT';  value5 = {zeros(48000,4)};
field6 = 'measurement_Iteration'; value6 = 1;

% Inititalizing BTE struct
BRIR_Data = struct(field1,value1,field2,value2,field3,value3,field4,...
    value4,field5,value5,field6,value6);

% Clearing unneeded initiazation variables
clear field1 field2 field3 field4 field5 field6 value1 value2 value3...
    value4 value5 value6 struct_Init

% Loading in BTE (measurement parameter) Data as a cell array
load(param_Str);

% For loop loading all of the raw BTE data into the struct
for index = 1:n_Meas
    
    try
        data_TEMP = load(strcat(data_Str,BRIR_Param_Data{index,4},'.mat'));
    catch
        disp(strcat("The file ",BRIR_Param_Data{index,4},'.mat',' does not exist'));
        continue
    end
    
    BRIR_Data(index).abscissa = BRIR_Param_Data{index,1};
    BRIR_Data(index).ordinate = BRIR_Param_Data{index,2};
    BRIR_Data(index).rotation = BRIR_Param_Data{index,3};
    BRIR_Data(index).IR_LEFT = data_TEMP.data(1).IR;
    BRIR_Data(index).IR_RIGHT = data_TEMP.data(2).IR;
end
