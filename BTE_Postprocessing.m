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

struct_Init = cell(1,num_Measurements);

for index=1:num_Measurements
   struct_Init{index} = 0;
end

% Creating the fields and initializing the values for the BTE Data Struct

field1 = 'abscissa';  value1 = struct_Init;
field2 = 'ordinate';  value2 = struct_Init;
field3 = 'rotation';  value3 = struct_Init;
field4 = 'IR';  value4 = {zeros(48000,4)};
field5 = 'measurement_Number'; value5 = 1;

BTE_Data = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);
clear field1 field2 field3 field4 field5 value1 value2 value3 value4 value5...
    struct_Init

% for index = 1 : k(1)
% variable(index) = load(sprintf('internalforce_%d',index))
% end