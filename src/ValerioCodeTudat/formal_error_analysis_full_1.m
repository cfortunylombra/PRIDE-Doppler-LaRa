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
% This script allows to calculate and plot the formal errors of the
% estimated parameters defined in the TUDAT file LaRa_POD.cpp as function of
% the tracking cadence, the number of stations and the observation
% correlation coefficient.
%
% 
% Inputs:
%    Use the otputs provided by the TUDAT file LaRa_POD.cpp (.dat files)
%    Add a .tex file with the names of the stations involved
%
%
% Other m-files required: none
% Subfunctions: multi_station_observation_weight.m
% MAT-files required: none
%

% Author: Valerio Filice
% email address: filicevalerio@gmail.com  
% Last revision: 20-Oct-2019

%% DEFINE DATA DIRECTORY AND LOAD FILES

format long; clear; close all; clc;

% Path of all the data
dataDirectory = strrep(matlab.desktop.editor.getActiveFilename,'src\ValerioCodeTudat\formal_error_analysis_full_1.m','output\POD_LaRa\');

% Load estimation information matrix
estimationInformationMatrix = load( strcat( dataDirectory, 'estimation_information_matrix.dat' ));

% Load information matrix transformation diagonal
informationMatrixTransformationDiagonal = load( strcat( dataDirectory,...
    'estimation_information_matrix_normalization.dat' ));

% Load overall observation time
observationTimes = load( strcat( dataDirectory, 'concatenated_times.dat' ));

% Load observation link ends
observationLinkEnds = load( strcat( dataDirectory, 'concatenated_link_ends.dat' ));

% List of Ground Station Names
groundStationID = ["BADARY"; "CEDUNA"; "HARTRAO"; "HART15M"; "HOBART12"; "HOBART26"; "TIANMA65"; "WARK30M"; "EFLSBERG"; "IRBENE"; "YEBES40M"; "MEDICINA"; "WETTZELL"; "ONSALA60"; "WRT0"];

%% Sort Data

% Define the transmitter antenna
index_transmitter = 0;

% 1 day in seconds
one_day = 86400; %seconds

% Define 1 mas
mas = pi / ( 180.0 * 1000.0 * 3600.0 );

% Extract the transmitter observation times
transmitterObservationTimes = observationTimes( observationLinkEnds == index_transmitter );

% Extract the transmitter information matrix
transmitterInformationMatrix = estimationInformationMatrix( observationLinkEnds == index_transmitter , :);

% Order the overall observation times, and subtract with the first time element
observationTimes( observationLinkEnds == index_transmitter ) = [];
observationTimes = [ transmitterObservationTimes; observationTimes];
observationTimes = (observationTimes - observationTimes(1) );

% Substract the transmitter observation times with the first time element
transmitterObservationTimes = transmitterObservationTimes - transmitterObservationTimes( 1 );

% Order the estimation information matrix
estimationInformationMatrix( observationLinkEnds == index_transmitter , :) = [];
estimationInformationMatrix = [ transmitterInformationMatrix ; estimationInformationMatrix ];

% Order the observation link-ends, and add +1 to all the elements
observationLinkEnds( observationLinkEnds == index_transmitter ) = [];
observationLinkEnds = [ones( [ length( transmitterObservationTimes ) 1]); 
    observationLinkEnds( observationLinkEnds > index_transmitter ) + 1];

% Add the name of the transmitter to the ground station names
groundStationID = [ 'DSS63' ; groundStationID];

%% Define Inverse Apriopri Covariance 

% Define the sigma apriori
sigmaAPriori = [ones(1,3)*1000, ones(1,3)*0.0002, 0.014, deg2rad(0.075)/one_day, ones(1,3)*30 ...
                23*mas, 26*mas, 22*mas, 22*mas, 18*mas, 19*mas, 16*mas, 16*mas, ones(1,20)*2*mas];

% Define the inverse apriori covariance from the sigma apriori            
inverseAPrioriCovariance = diag(1./sigmaAPriori.^2);

% Normalize the inverse apriori covariance matrix
normalizedInverseAprioriCovarianceMatrix = zeros(size(inverseAPrioriCovariance));
for index=1:length(informationMatrixTransformationDiagonal)
    for ind=1:length(informationMatrixTransformationDiagonal)
        normalizedInverseAprioriCovarianceMatrix(index,ind) = ...
            inverseAPrioriCovariance(index,ind) / (informationMatrixTransformationDiagonal(index)*...
            informationMatrixTransformationDiagonal(ind));
    end
end

%% DEFINE VARIANCE-COVARIANCE MATRIX PARAMETERS

% Speed of light [m/s]
c = 299792458.0;

% Doppler noise 
sigma_PRIDE = 0.05E-3 / c ;
sigma_DSS63 = 0.05E-3 / c ;

% Define variance
variance = [ sigma_DSS63 repelem( sigma_PRIDE, length(groundStationID)-1)];

% Define possible correlation coefficients
correlationCoefficient = (0:0.1:0.9);

%% FORMAL ERRORS

% Take the observations from the non-transmitter antennas
prideObservationTimes = observationTimes( observationLinkEnds > 1 );
prideObservationLinkEnds = observationLinkEnds( observationLinkEnds > 1 );
prideInformationMatrix = estimationInformationMatrix( observationLinkEnds > 1 , : );

weekPerMonth = 4; %Select the number of weeks per month 
for week = 1 : weekPerMonth 
    fprintf('Week:')
    disp(week)
        
    % Create the observation times
    obsTime = [];
    for numberOfMonths = 1 : 25
        for numberOfWeek = 1 : week
            for numberPerWeek = 1 : 2
                for ind = 1 : 2   
                    if ind == 1
                        obsTime = [ obsTime, (numberOfMonths - 1) * 28 * one_day + (numberOfWeek - 1) * 7 * one_day + ...
                            (numberPerWeek - 1) * 3.25 * one_day - one_day ];
                    else
                        obsTime = [ obsTime, (numberOfMonths - 1) * 28 * one_day + (numberOfWeek - 1) * 7 * one_day + ...
                            (numberPerWeek - 1) * 3.25 * one_day + 2*one_day ];
                    end                    
                end
            end
        end
    end
    
    % Create the selected observation times and link-ends
    prideSelectedObservationTimes = [];
    prideSelectedObservationLinkEnds = [];
    prideSelectedInformationMatrix = [];
    
    % Iterate to each non-transmitter antenna
    for groundStationIDindex = 2 : length( groundStationID )
        currentStationTimes = prideObservationTimes( prideObservationLinkEnds == groundStationIDindex );
        currentStationLinkEnds = prideObservationLinkEnds( prideObservationLinkEnds == groundStationIDindex );
        currentStationInformationMatrix = prideInformationMatrix( prideObservationLinkEnds == groundStationIDindex, : );
        
        % Select the observation times, link ends, information matrix inside the time spans
        for index = 2 : 2 : length( obsTime )
            prideSelectedObservationTimes = [ prideSelectedObservationTimes;  ...
                currentStationTimes( currentStationTimes >= obsTime(index - 1) & currentStationTimes <= obsTime( index ))];
            prideSelectedObservationLinkEnds = [ prideSelectedObservationLinkEnds; ...
                currentStationLinkEnds( currentStationTimes >= obsTime(index - 1) & currentStationTimes <= obsTime( index ))];
            prideSelectedInformationMatrix = [ prideSelectedInformationMatrix; ...
                currentStationInformationMatrix( currentStationTimes >= obsTime(index - 1) & currentStationTimes <= obsTime( index ),:)];
        end
    end
    
    % Add the information from the transmitter antenna
    time = [transmitterObservationTimes; prideSelectedObservationTimes];
    informationMatrix = [ transmitterInformationMatrix; prideSelectedInformationMatrix];
    linkEnds = [ones( [ length( transmitterObservationTimes ) 1]); prideSelectedObservationLinkEnds];
    
    % Create observation times, link-ends and estimation information matrix
    obsTimes = [];
    obsLinkEnds = [];
    estInformationMatrix = [];
    
    % Iterate for each ground station
    for groundStationIDindex = 1 : length( groundStationID )
        fprintf('groundStationIDindex:')
        disp(groundStationIDindex)
        
        % Append the obsTimes, obslinkEnds and estInformationMatrix for each ground station
        obsTimes = [obsTimes; time( linkEnds == groundStationIDindex )];
        obsLinkEnds = [obsLinkEnds; linkEnds( linkEnds == groundStationIDindex )];
        estInformationMatrix = [estInformationMatrix; informationMatrix( linkEnds == groundStationIDindex, : )];
        
        % Iterate for each correlation coefficient
        for correlationCoefficientIndex = 1 : length(correlationCoefficient)
            fprintf('correlationCoefficientIndex:')
            disp(correlationCoefficientIndex)    
            
            % Set weights
            [ H, Wgt ] = multi_station_observation_weight( obsTimes, obsLinkEnds, estInformationMatrix,...
                variance, correlationCoefficient( correlationCoefficientIndex ));
            
            % Perform formal error analysis
            inverseNormalizedCovarianceMatrix = H'*Wgt*H + normalizedInverseAprioriCovarianceMatrix;
            
            % Invert
            normalizedCovarianceMatrix = inv(inverseNormalizedCovarianceMatrix);
            
            % Create the unnormalized covariance matrix
            unnormalizedCovarianceMatrix = zeros(size(normalizedCovarianceMatrix));
            for index=1:length(informationMatrixTransformationDiagonal)
                for ind=1:length(informationMatrixTransformationDiagonal)
                    unnormalizedCovarianceMatrix(index,ind) = ...
                        normalizedCovarianceMatrix(index,ind) / (informationMatrixTransformationDiagonal(index) * ...
                        informationMatrixTransformationDiagonal(ind));
                end
            end
            
            % Take the standard variation
            sigma0 = sqrt(diag(unnormalizedCovarianceMatrix));
            
            sigmas_mars_state_comb(:, correlationCoefficientIndex, groundStationIDindex)  = sigma0(1 : 6);
            
            sigmas_comb(:, correlationCoefficientIndex, groundStationIDindex)  = sigma0(7 : end);
        end
    end
    sigmas_comb_cell{ week } = sigmas_comb;
end

%% Plot

load('formal1.mat')

col = hsv(length( groundStationID ));

% The aim is to understand whether the sigmas are reduced when the amount of data is maximum
figure( 1 )
for week = 1 : 4
   subplot(2, 2, week)
   hold on
   current_sigmas_comb = sigmas_comb_cell{ week };
   for gsIndex = 1 : length(groundStationID)
       for corrIndex = 1 : length( correlationCoefficient )
           plot(gsIndex, current_sigmas_comb( 1, corrIndex, gsIndex ), ':s', 'Color', col( corrIndex, : ))
       end
   end
   xlabel('Number of stations')
   ylabel('$\sigma_{F}$','Interpreter','latex','FontSize',19)
   string = sprintf('Tracking sessions: %u week per month', week);
   title(string)
   box on
   grid on
end
legend('\rho = 0','\rho = 0.1','\rho = 0.2','\rho = 0.3','\rho = 0.4','\rho = 0.5','\rho = 0.6','\rho = 0.7',...
    '\rho = 0.8','\rho = 0.9', 'Position',[0.926953125 0.42755680598996 0.05390625 0.168323863636364])
saveas(gcf,'sigma_weeks.png')