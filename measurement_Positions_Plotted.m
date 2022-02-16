%% measurement_Positions_Plotted.m
%
% This script works through﻿Kabzinski's Time-LS CTC results
%
% Last updated by Michael Klunk 02/02/2021

%% Initializing the Script

clear
close all
clc
load('BTE_Locations.mat');

%% Positions of the Loudspeakers

% Left speaker (facing the TV)
speaker_Pos(1,:) = [0.561, (3.475-0.576)];

% Right speaker (facing the TV)
speaker_Pos(2,:) = [(2.930-0.559), (3.475-0.556)];

%% Plotting the Postions of Measurements

% Creating the figure
figure

    % MEASUREMENT POSITIONS
    scatter(BTE_Locations(:,6), BTE_Locations(:,7), 150, 'd', 'filled');
    grid on, grid minor
    axis equal
    xlim([0,2.93]), ylim([0,3.475])
    
    % Speaker Positions
    hold on 
    scatter(speaker_Pos(:,1), speaker_Pos(:,2), 250, 'vg', 'LineWidth', 1.5);
    xlabel('Width of Audio-Booth (m)')
    ylabel('Length of Audio-Booth (m)')
    title('Loudspeaker and H.A.T.S. Measurement Positions')



