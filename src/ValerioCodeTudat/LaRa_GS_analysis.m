%% Copyright (c) 2022, Carlos Fortuny Lombraña - Valerio Filice - Dominic Dirkx
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%% Description
% This script allows to calculate and plot the mean elevation angle of the
% ground stations and of the lander. Moreover, it plots the number of
% possible observations of the receivers.
%
% 
% Inputs:
%    Use the otputs provided by the TUDAT file LaRa_GS_analysis.cpp (.dat files)
%    Add a .tex file with the names of the stations involved
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Carlos Fortuny Lombraña
% email address: C.FortunyLombrana@student.tudelft.nl  

%% DEFINE DATA DIRECTORY AND LOAD FILES

format long; clear; close all; clc;

% Path of all the data
dataDirectory = strrep(matlab.desktop.editor.getActiveFilename,'src\ValerioCodeTudat\LaRa_GS_analysis.m','output\GS_LaRa\');

% List of Ground Station Names
groundStationNames = [ "DSS63"; "BADARY"; "CEDUNA"; "HARTRAO"; "HART15M"; "HOBART12"; "HOBART26"; "TIANMA65"; "WARK30M"; "EFLSBERG"; "IRBENE"; "YEBES40M"; "MEDICINA"; "WETTZELL"; "ONSALA60"; "WRT0"];

% Load Overall Observation Time
EphemerisTime = load(strcat(dataDirectory,'observation_time.dat'));
start_epoch = EphemerisTime(1);
EphemerisTime = EphemerisTime - EphemerisTime( 1 );

% Load DSS63 Observation Time
DSS63ObservationTime = load(strcat(dataDirectory,'DSS63_observation_time.dat'));
DSS63ObservationTime = DSS63ObservationTime - DSS63ObservationTime( 1 );

% Load Ground Observation times and complement with DSS63 Observation Time
groundStationObservationTimes = load(strcat(dataDirectory,'ground_station_observation_time.dat'));
groundStationObservationTimes = groundStationObservationTimes - groundStationObservationTimes( 1 );
groundStationObservationTimes = [ DSS63ObservationTime; groundStationObservationTimes ];

% Load DSS63 Elevation Angles [rad]
DSS63Elevation = load(strcat(dataDirectory,'DSS63_elevation.dat'));

% Load Earth Elevation Angles [rad]
earthElevation = load(strcat(dataDirectory,'earth_elevation.dat'));

% Load Earth Azimuth Angles [rad]
earthAzimuth = load(strcat(dataDirectory,'earth_azimuth.dat'));

% Load Ground Station Elevations and complement with DSS63 Observation Time
groundStationElevations = load(strcat(dataDirectory,'ground_station_elevation.dat'));
groundStationElevations = [ DSS63Elevation; groundStationElevations ];

% Load the Ground Station IDs and complement with DSS63
groundStationIDs = load(strcat(dataDirectory,'ground_station_ids.dat')) + 1;
groundStationIDs = [zeros([length( DSS63Elevation ) 1]); groundStationIDs ] + 1;

% 1 day in seconds
one_day = 86400; %seconds

% For plotting, just show the observations each week
step_days = 7;

% Viability setting
max_elevation = deg2rad(20);

% J2000 epoch
J2000_in_Julian_days = 2451545.0;

% Save date when the simulation starts
start_date = datetime(start_epoch/one_day+J2000_in_Julian_days, 'convertfrom','juliandate');
%% Elevation and Azimuth Angles of the Earth as seen by LaRa

figure(1)
xgscatter = rad2deg(earthAzimuth(earthElevation >= 0));
ygscatter = rad2deg(earthElevation(earthElevation >= 0));
groupgscatter = ygscatter >= 35 & ygscatter <= 45;
gscatter(xgscatter, ygscatter, groupgscatter, 'br', 'o', 5, 'off')
axis([-180 180 0 90])
set(gca,'XTick',(-180 : 45 : 180));
xlabel('Azimuth angle [deg]')
ylabel('Elevation angle [deg]')
grid on
sgtitle(join(['Mission Duration: 700 days - Start Date: ',datestr(start_date)]))
legend('No Tracking Windows of LaRa Found','Tracking Windows of LaRa','Location','northeastoutside')
set(gcf,'Position', 0.85*get(0, 'Screensize'));
set(gca,'FontSize',14);

%% Transmitter Elevation History, Mean Elevation

figure(2)
set(gca,'FontSize',20)
hold on

% Polynomial curve fitting (order 4) 
[fit, S, mu] = polyfit(DSS63ObservationTime, DSS63Elevation, 4);

% Polyval is the polynomial evaluation
plot((0:(step_days*one_day):EphemerisTime(end)) / one_day , rad2deg( polyval( fit,(0:(step_days*one_day):EphemerisTime(end)) , [], mu )),'DisplayName','Mean Elevation Angle')
plot( DSS63ObservationTime / one_day, rad2deg( DSS63Elevation), 'o','DisplayName','Elevation Angle')

hold off

xlabel('Mission Time [days]')
ylabel('Elevation Angle [deg]')
box on
grid on
sgtitle(join(['Start Date: ',datestr(start_date)]))
legend('Location','northeastoutside')
set(gcf,'Position', 0.85*get(0, 'Screensize'));
set(gca,'FontSize',14);

%% Receiver Mean Elevation

figure(3)
set(gca,'FontSize',20)
set(0,'defaultaxeslinestyleorder',{'-*','-+','-o','-x','-s'})
hold on
groundStationNames_legend = [];
for groundStationIDindex = 1 : length( groundStationNames )
    
    % Angle viability due to Sun
    % If it is not satisfied, the next ground station is evaluated
    if isempty(groundStationObservationTimes( groundStationIDs == groundStationIDindex & groundStationElevations >= max_elevation))
        continue
    end
    
    % Plot only the ground stations visible
    groundStationNames_legend = [groundStationNames_legend; groundStationNames(groundStationIDindex)];
    
    % Polynomial curve fitting (order 4) 
    [fit, S, mu] = polyfit(groundStationObservationTimes( groundStationIDs == groundStationIDindex ), groundStationElevations( groundStationIDs == groundStationIDindex ), 4);
    
    plot((0:(step_days*one_day):EphemerisTime(end))/one_day, rad2deg( polyval( fit,0:(step_days*one_day):EphemerisTime(end), [], mu )))    
end
hold off
ylim([0 90])
xlabel('Mission time [days]')
ylabel('Elevation angle [deg]')
sgtitle(join(['Start Date: ',datestr(start_date)]))
legend( groundStationNames_legend, 'Location', 'northeastoutside')
box on
grid on
set(gcf,'Position', 0.85*get(0, 'Screensize'));

%% Number of Possible Observations for each ground station
        
figure(4)
subplot(1,2,1)
hold on
groundStationNumberObservations = [];

for groundStationIDindex = 1 : length(groundStationNames)
    
    %Check the viability setting 
    currentGroundStationObservationTimes = groundStationObservationTimes( groundStationIDs == groundStationIDindex & groundStationElevations >= max_elevation);
    
    % Counting the number of observables using the length function
    groundStationNumberObservations = [ groundStationNumberObservations length(currentGroundStationObservationTimes)];
    
    % Plot
    scatter( currentGroundStationObservationTimes/one_day, repmat( groundStationIDindex, 1, length(currentGroundStationObservationTimes)), 200, '.')
end
hold off
ay = gca;
ay.YTick = 1:numel(groundStationNames);
ay.YTickLabel = groundStationNames;
ay.YLim = [0 numel(groundStationNames)+1];
ay.XLim = [-10 inf];
ay.XLabel.String = 'Mission time [days]';
ay.XMinorGrid = 'on';

subplot(1,2,2)
barh( categorical( groundStationNames, groundStationNames ), groundStationNumberObservations )
text( groundStationNumberObservations + 5 , 1:length( groundStationNumberObservations ), ...
    num2str( groundStationNumberObservations' ),'vert','middle','horiz','left');
box off
ax = gca;
xlim([ -5 ax.XLim( 2 ) + 500])
xlabel('Number of Observations')
sgtitle(join(['Start Date: ',datestr(start_date)]))
set(gcf,'Position', 0.85*get(0, 'Screensize'));
