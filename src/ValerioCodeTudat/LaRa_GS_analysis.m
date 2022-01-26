%% Copyright (c) 2019, Valerio Filice - Dominic Dirkx
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

% Author: Valerio Filice
% email address: filicevalerio@gmail.com  
% Last revision: 20-Oct-2019

%% DEFINE DATA DIRECTORY AND LOAD FILES

format long; clear; close all; clc;

% Path of all the data
dataDirectory = strrep(matlab.desktop.editor.getActiveFilename,'src\ValerioCodeTudat\LaRa_GS_analysis.m','output\GS\');

% List of Ground Station Names
groundStationNames = [ "DSS63"; "BADARY"; "CEDUNA"; "HARTRAO"; "HART15M"; "HOBART12"; "HOBART26"; "TIANMA65"; "WARK30M"; "EFLSBERG"; "IRBENE"; "YEBES40M"; "MEDICINA"; "WETTZELL"; "ONSALA60"; "WRT0"];

% Load Overall Observation Time
EphemerisTime = load(strcat(dataDirectory,'observation_time.dat'));
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

%% Elevation and Azimuth Angles of the Earth as seen by LaRa

figure(1)
%xgscatter = rad2deg(wrapToPi(earthAzimuth(earthElevation >= 0) + pi/2 ));
xgscatter = rad2deg(earthAzimuth(earthElevation >= 0));
ygscatter = rad2deg(earthElevation(earthElevation >= 0));
groupgscatter = ygscatter >= 35 & ygscatter <= 45;
gscatter(xgscatter, ygscatter, groupgscatter, 'br', 'o', 5, 'off')
axis([-180 180 0 90])
set(gca,'XTick',(-180 : 45 : 180));
xlabel('Azimuth angle [deg]')
ylabel('Elevation angle [deg]')
set(gcf,'Position', get(0, 'Screensize'));
        
%% Transmitter Elevation History

figure(2)
set(gca,'FontSize',20)
hold on

% Polynomial curve fitting 
[fit, S, mu] = polyfit(DSS63ObservationTime, DSS63Elevation, 4);

% 86400 sec = 1 day
% Plotting with steps of 12.5
% Polyval is the polynomial evaluation
plot((0:(12.5*86400):EphemerisTime(end)) / 86400 , rad2deg( polyval( fit,(0:(12.5*86400):EphemerisTime(end)) , [], mu )))
plot( DSS63ObservationTime / 86400, rad2deg( DSS63Elevation), 'o')
% plot( DSS63ObservationTime( DSS63Elevation >= 0) / 86400, rad2deg( DSS63Elevation( DSS63Elevation >= 0)), 'o')

hold off

xlabel('Mission Time [days]')
ylabel('Elevation Angle [deg]')
box on
set(gcf,'Position', get(0, 'Screensize'));

% Mean Elevation angle

figure( 3 )
set(gca,'FontSize',20)
% col = jet(length( groundStationID ));
set(0,'defaultaxeslinestyleorder',{'-*','-+','-o','-x','-s'})
hold on
for groundStationIDindex = 2 : length( groundStationNames )
    
    % Angle viability due to Sun
    % If it is not satisfied, the next ground station is evaluated
    if isempty(groundStationObservationTimes( groundStationIDs == groundStationIDindex & groundStationElevations >= deg2rad( 20 )))
        
        continue
        
    end
    
%     [fit, S, mu] = polyfit(groundStationObservationTimes( groundStationIDs == groundStationIDindex ...
%         & groundStationElevations >= deg2rad( 20 )), ...
%         groundStationElevations( groundStationIDs == groundStationIDindex & groundStationElevations >= deg2rad( 20 )), 4);
%
    [fit, S, mu] = polyfit(groundStationObservationTimes( groundStationIDs == groundStationIDindex ), ...
        groundStationElevations( groundStationIDs == groundStationIDindex ), 4);
    
    plot((0:(12.5*86400):EphemerisTime(end)) / 86400, rad2deg( polyval( fit,0:(12.5*86400):EphemerisTime(end), [], mu )))    
    
end
hold off
ylim([0 90])
xlabel('Mission time [days]')
ylabel('Elevation angle [deg]')
legend( groundStationNames(2:end) , 'Location', 'best')
box on
set(gcf, 'Position', get(0, 'Screensize'));

%% Number of Possible Observations
        
for Set = 16 : 16 : length( groundStationNames )
        
    if Set == 16
        
        groundStationID = {groundStationNames{ 1 : Set }};
        startIndex = 1;
        endIndex = 16;
        
    else
        
        groundStationID = {groundStationNames{ Set - 15 : Set }};
        startIndex = endIndex + 1 ;
        endIndex = startIndex + 15 ;
        
    end
    
    figure
    subplot(1,2,1)
    hold on
    groundSationNumberObservations = [];

    countFor = 1;
    for groundStationIDindex = startIndex : endIndex

        currentGroundStationObservationTimes = ...
            groundStationObservationTimes( groundStationIDs == groundStationIDindex & groundStationElevations >= deg2rad( 0 ));
        groundSationNumberObservations = [ groundSationNumberObservations length( currentGroundStationObservationTimes )];
        
        scatter( currentGroundStationObservationTimes / 86400, ...
            repmat( countFor, 1, length( currentGroundStationObservationTimes )), 'o')
        countFor = countFor + 1;
    end
    hold off
    ay = gca;
    ay.YTick = 1:numel(groundStationID);
    ay.YTickLabel = groundStationID;
    ay.YLim = [0 numel(groundStationID)+1];
    ay.XLim = [-10 inf];
    ay.XLabel.String = 'Observation Time [gg]';
    % ay.XGrid = 'on';
    ay.XMinorGrid = 'on';
    
    subplot(1,2,2)
    barh( categorical( groundStationID, groundStationID ), groundSationNumberObservations )
    text( groundSationNumberObservations + 5 , 1:length( groundSationNumberObservations ), ...
        num2str( groundSationNumberObservations' ),'vert','middle','horiz','left');
    box off
    ax = gca;
    xlim([ -5 ax.XLim( 2 ) + 500])
    xlabel('Number of Observations')
    set(gcf, 'Position', get(0, 'Screensize'));
        
end

if Set < length( groundStationNames )
            
    groundStationID = {groundStationNames{ Set : end }};
    startIndex = Set ;
    endIndex = length( groundStationNames ) ;
    
        figure
    subplot(1,2,1)
    hold on
    groundSationNumberObservations = [];

    countFor = 1;
    for groundStationIDindex = startIndex : endIndex

        currentGroundStationObservationTimes = ...
            groundStationObservationTimes( groundStationIDs == groundStationIDindex & groundStationElevations >= deg2rad( 20 ));
        groundSationNumberObservations = [ groundSationNumberObservations length( currentGroundStationObservationTimes )];
        
        scatter( currentGroundStationObservationTimes / 86400, ...
            repmat( countFor, 1, length( currentGroundStationObservationTimes )), 'o')
        countFor = countFor + 1;
    end
    hold off
    ay = gca;
    ay.YTick = 1:numel(groundStationID);
    ay.YTickLabel = groundStationID;
    ay.YLim = [0 numel(groundStationID)+1];
    ay.XLim = [-10 inf];
    ay.XLabel.String = 'Observation Time [gg]';
    % ay.XGrid = 'on';
    ay.XMinorGrid = 'on';
    
    subplot(1,2,2)
    barh( categorical( groundStationID, groundStationID ), groundSationNumberObservations )
    text( groundSationNumberObservations + 5 , 1:length( groundSationNumberObservations ), ...
        num2str( groundSationNumberObservations' ),'vert','middle','horiz','left');
    box off
    ax = gca;
    xlim([ -50 ax.XLim( 2 ) + 500])
    xlabel('Number of Observations')
    set(gcf, 'Position', get(0, 'Screensize'));
    
end

    
