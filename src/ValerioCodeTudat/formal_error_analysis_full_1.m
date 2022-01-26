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

dataDirectory = '/Users/valeriofilice/Tudat/tudatBundle/tudatApplications/Outputs/LaRa_POD_full_10gs/';

estimationInformationMatrix = load( strcat( dataDirectory, 'EstimationInformationMatrix.dat' ));

informationMatrixTransformationDiagonal = load( strcat( dataDirectory,...
    'EstimationInformationMatrixNormalization.dat' ));

observationTimes = load( strcat( dataDirectory, 'ObservationTimes.dat' ));

observationLinkEnds = load( strcat( dataDirectory, 'ObservationLinkEnds.dat' ));

groundStationID = textscan( fopen( strcat(dataDirectory, 'gs_location.txt') ), '%s','delimiter','\n');

%% Sort Data

transmitterObservationTimes = observationTimes( observationLinkEnds == 3 );

transmitterInformationMatrix = estimationInformationMatrix( observationLinkEnds == 3 , :);

observationTimes( observationLinkEnds == 3 ) = [];
observationTimes = [ transmitterObservationTimes; observationTimes];
observationTimes = (observationTimes - observationTimes(1) );

transmitterObservationTimes = transmitterObservationTimes - transmitterObservationTimes( 1 );

estimationInformationMatrix( observationLinkEnds == 3 , :) = [];
estimationInformationMatrix = [ transmitterInformationMatrix ; estimationInformationMatrix ];

observationLinkEnds( observationLinkEnds == 3 ) = [];
observationLinkEnds = [ones( [ length( transmitterObservationTimes ) 1]); observationLinkEnds( observationLinkEnds < 3 ) + 2; ...
    observationLinkEnds( observationLinkEnds > 3 ) + 1 ];

groundStationID = [ 'DSS63' ; groundStationID{ 1 }];

milliArcSecondToRadian = pi / ( 180.0 * 1000.0 * 3600.0 );

%% Define Inverse Apriopri Covariance

FCN_0 = 2*pi/(-240*86400); %[rad/dd]
sigma_FCN_0 = 2 * pi /(240*86400)^2 * 50 * 86400; %[rad/s]

sigmaAPriori = [0.07 sigma_FCN_0 10 10 10000 repelem(milliArcSecondToRadian * 15, 8) repelem(10 * milliArcSecondToRadian, 8) ...
    repelem(50*milliArcSecondToRadian, 4)]';

inverseAPrioriCovariance = ...
    diag([1 1 1 1 1 1 ...
    1/(0.07)^2 1/(sigma_FCN_0)^2 1/(10)^2 1/(10)^2 1/(10000)^2 ...
    1/(milliArcSecondToRadian * 15)^2 1/(milliArcSecondToRadian * 15)^2 1/(milliArcSecondToRadian * 15)^2 ...
    1/(milliArcSecondToRadian * 15)^2 1/(milliArcSecondToRadian * 15)^2 1/(milliArcSecondToRadian * 15)^2 ...
    1/(milliArcSecondToRadian * 15)^2 1/(milliArcSecondToRadian * 15)^2 ...
    1/(10 * milliArcSecondToRadian )^2 1/(10 * milliArcSecondToRadian )^2 ...
    1/(10 * milliArcSecondToRadian )^2 1/(10 * milliArcSecondToRadian )^2 ...
    1/(10 * milliArcSecondToRadian )^2 1/(10 * milliArcSecondToRadian )^2 ...
    1/(10 * milliArcSecondToRadian )^2 1/(10 * milliArcSecondToRadian )^2 ...
    1/(50*milliArcSecondToRadian)^2 1/(50*milliArcSecondToRadian)^2 ...
    1/(50*milliArcSecondToRadian)^2 1/(50*milliArcSecondToRadian)^2]);

% inverseAPrioriCovariance = zeros( length( informationMatrixTransformationDiagonal ));

normalizedInverseAprioriCovarianceMatrix = zeros(size(inverseAPrioriCovariance));
for index=1:length(informationMatrixTransformationDiagonal)
    
    for ind=1:length(informationMatrixTransformationDiagonal)
        
        normalizedInverseAprioriCovarianceMatrix(index,ind) = ...
            inverseAPrioriCovariance(index,ind) / (informationMatrixTransformationDiagonal(index) * ...
            informationMatrixTransformationDiagonal(ind));
    end
    
end

%% DEFINE VARIANCE-COVARIANCE MATRIX PARAMETERS

c = 299792458.0;
sigma_PRIDE = 0.05E-3 / c ;
sigma_DSS63 = 0.05E-3 / c ;

variance = [ sigma_DSS63 repelem( sigma_PRIDE, 9 )];

correlationCoefficient = (0:0.1:0.9);

%% FORMAL ERRORS

%Select frequency of the observations for PRIDE

prideObservationTimes = observationTimes( observationLinkEnds > 1 );
prideObservationLinkEnds = observationLinkEnds( observationLinkEnds > 1 );
prideInformationMatrix = estimationInformationMatrix( observationLinkEnds > 1 , : );

matrixArrayIndex = 1;
weekPerMonth = 4; %Select the number of weeks per month 
for week = 1 : weekPerMonth 
    obsTime = [];
    for numberOfMonths = 1 : 25
        for numberOfWeek = 1 : week
            for numberPerWeek = 1 : 2
                for ind = 1 : 2
                    
                    if ind == 1
                        obsTime = [ obsTime (numberOfMonths - 1) * 28 * 86400 + (numberOfWeek - 1) * 7 * 86400 + ...
                            (numberPerWeek - 1) * 3.25 * 86400 - 86400 ];
                    else
                        obsTime = [ obsTime (numberOfMonths - 1) * 28 * 86400 + (numberOfWeek - 1) * 7 * 86400 + ...
                            (numberPerWeek - 1) * 3.25 * 86400 + 2*86400 ];
                    end
                    
                end
            end
        end
    end
    
    prideSelectedObservationTimes = [];
    prideSelectedObservationLinkEnds = [];
    prideSelectedInformationMatrix = [];
    for groundStationIDindex = 2 : length( groundStationID )
        
        currentStationTimes = prideObservationTimes( prideObservationLinkEnds == groundStationIDindex );
        currentStationLinkEnds = prideObservationLinkEnds( prideObservationLinkEnds == groundStationIDindex );
        currentStationInformationMatrix = prideInformationMatrix( prideObservationLinkEnds == groundStationIDindex, : );
        
        for index = 2 : 2 : length( obsTime )
            
            prideSelectedObservationTimes = [ prideSelectedObservationTimes;  ...
                currentStationTimes( currentStationTimes >= obsTime(index - 1) & currentStationTimes <= obsTime( index ))];
            prideSelectedObservationLinkEnds = [ prideSelectedObservationLinkEnds; ...
                currentStationLinkEnds( currentStationTimes >= obsTime(index - 1) & currentStationTimes <= obsTime( index ))];
            prideSelectedInformationMatrix = [ prideSelectedInformationMatrix; ...
                currentStationInformationMatrix( currentStationTimes >= obsTime(index - 1) & currentStationTimes <= obsTime( index ),:)];
        end
    end
    
    time = [transmitterObservationTimes; prideSelectedObservationTimes];
    informationMatrix = [ transmitterInformationMatrix; prideSelectedInformationMatrix];
    linkEnds = [ones( [ length( transmitterObservationTimes ) 1]); prideSelectedObservationLinkEnds];
    
    obsTimes = [];
    obsLinkEnds = [];
    estInformationMatrix = [];
    
    for groundStationIDindex = 1 : length( groundStationID )
        
        obsTimes = [obsTimes; time( linkEnds == groundStationIDindex )];
        
        obsLinkEnds = [obsLinkEnds; linkEnds( linkEnds == groundStationIDindex )];
        
        estInformationMatrix = [estInformationMatrix; informationMatrix( linkEnds == groundStationIDindex, : )];
        
        for correlationCoefficientIndex = 1 : length(correlationCoefficient)
            
            [ H, Wgt ] = multi_station_observation_weight( obsTimes, obsLinkEnds, estInformationMatrix,...
                variance, correlationCoefficient( correlationCoefficientIndex ));
            
            % PERFORM FORMAL ERROR ANALYSIS
            
            inverseNormalizedCovarianceMatrix = H'*Wgt*H + normalizedInverseAprioriCovarianceMatrix;
            
            normalizedCovarianceMatrix = inv(inverseNormalizedCovarianceMatrix);
            
            unnormalizedCovarianceMatrix = zeros(size(normalizedCovarianceMatrix));
            for index=1:length(informationMatrixTransformationDiagonal)
                
                for ind=1:length(informationMatrixTransformationDiagonal)
                    
                    unnormalizedCovarianceMatrix(index,ind) = ...
                        normalizedCovarianceMatrix(index,ind) / (informationMatrixTransformationDiagonal(index) * ...
                        informationMatrixTransformationDiagonal(ind));
                end
                
            end
            
            sigma0 = sqrt(diag(unnormalizedCovarianceMatrix));
            
            sigmas_mars_state_comb(:, correlationCoefficientIndex, groundStationIDindex)  = sigma0(1 : 6);
            
            sigmas_comb(:, correlationCoefficientIndex, groundStationIDindex)  = sigma0(7 : end);
                        
        end
        
    end
    
    sigmas_comb_cell{ week } = sigmas_comb;
    
end

%% Plot

col = hsv(length( groundStationID ));

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
   ylim([0.015 0.04])
   string = sprintf('Tracking sessions: %u week per month', week);
   title(string)
   box on
   grid on
       
end

legend('\rho = 0','\rho = 0.1','\rho = 0.2','\rho = 0.3','\rho = 0.4','\rho = 0.5','\rho = 0.6','\rho = 0.7',...
    '\rho = 0.8','\rho = 0.9', 'Position',[0.926953125 0.42755680598996 0.05390625 0.168323863636364])







