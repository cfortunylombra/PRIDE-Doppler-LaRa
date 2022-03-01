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
% the tracking cadence, the mission duration and the tracking technique. The 
% ratio between the formal error and the a priori uncertainties is plotted
% too.
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
dataDirectory = strrep(matlab.desktop.editor.getActiveFilename,'src\ValerioCodeTudat\formal_error_analysis_full_3.m','output\POD_LaRa\');

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
    observationLinkEnds( observationLinkEnds > index_transmitter ) + 1 ];

groundStationID = [ 'DSS63' ; groundStationID];

%% Inverse Apriopri Covariance

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
            inverseAPrioriCovariance(index,ind) / (informationMatrixTransformationDiagonal(index) * ...
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
variance = [ sigma_DSS63 repelem( sigma_PRIDE, length(groundStationID)-1 )];

% Define possible correlation coefficients
correlationCoefficient = (0:0.1:0.9);

%% PLOT FORMAL ERRORS - TRANSMITTER OBSERVATIONS ONLY / PRIDE OBSERVATIONS ONLY / TRANSMITTER + PRIDE OBSERVATIONS

% Obtain a large step for the analysis
step = divisors( length( transmitterObservationTimes ));
step = step( end - 2 );

% Array with the evaluation times
evaluationTimes_DSS63 = [];

% Create a transmitter weight matrix
transmitterWeightMatrix = eye( length( transmitterObservationTimes )) * 1 / sigma_DSS63^2;

% Iterate along the transmitter observation times with a step
indexIter = 1;
for Set = 0 : step : length( transmitterObservationTimes )
    
    % End condition
    if Set == length( transmitterObservationTimes )     
        currentObservationTime = transmitterObservationTimes( end );
        currentInformationMatrix = transmitterInformationMatrix; 
        currentWeightMatrix = transmitterWeightMatrix;
   
    % Initial condition
    elseif Set == 0   
        currentObservationTime = 0;
        currentInformationMatrix = 0;
        currentWeightMatrix = 0; 
        
    % Condition in-between
    else
        currentObservationTime = transmitterObservationTimes( Set );
        currentInformationMatrix = transmitterInformationMatrix( 1 : (Set + step - 1), : );
        currentWeightMatrix = transmitterWeightMatrix( 1 : (Set + step - 1), 1 : (Set + step - 1) ); 
    end
    
    % Create inverse normalized covariance matrix
    inverseNormalizedCovarianceMatrix = currentInformationMatrix' * currentWeightMatrix * currentInformationMatrix + ...
        normalizedInverseAprioriCovarianceMatrix;
    
    % Create normalized covariance matrix (Invert)
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
    
    % Take the standard deviations
    sigma0 = sqrt(diag(unnormalizedCovarianceMatrix))';
    
    sigmas_DSS63(:, indexIter)  = sigma0(7:end);
       
    evaluationTimes_DSS63 = [ evaluationTimes_DSS63 currentObservationTime ];

    indexIter = indexIter + 1;
end

% Take the observations from the non-transmitter antennas
prideObservationTimes = observationTimes( observationLinkEnds > 1 );
prideObservationLinkEnds = observationLinkEnds( observationLinkEnds > 1 );
prideInformationMatrix = estimationInformationMatrix( observationLinkEnds > 1 , : );

matrixArrayIndex = 1;
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
                        obsTime = [ obsTime (numberOfMonths - 1) * 28 * one_day + (numberOfWeek - 1) * 7 * one_day + ...
                            (numberPerWeek - 1) * 3.25 * one_day - one_day ];
                    else
                        obsTime = [ obsTime (numberOfMonths - 1) * 28 * one_day + (numberOfWeek - 1) * 7 * one_day + ...
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
        fprintf('groundStationIDindex:')
        disp(groundStationIDindex)
        
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
    
    % Iterate for each evaluated observation of the transmitter
    indexIter = 1;
    for Set = 1 : length( evaluationTimes_DSS63 )    
        fprintf('Set:')
        disp(Set)
        
        % End condition
        if Set == length( evaluationTimes_DSS63 )
            currentObservationTime = time;
            currentInformationMatrix = informationMatrix;
            currentObservationLinkEnds =  linkEnds;
        
        % Else condition
        else
            currentObservationTime = time( time <= evaluationTimes_DSS63( Set ));
            currentInformationMatrix = informationMatrix( time <= evaluationTimes_DSS63( Set ), : );
            currentObservationLinkEnds = linkEnds( time <= evaluationTimes_DSS63( Set ));
        end
        
        % Correlation coefficient chosen is 0.5, and the weighted matrices are created
        [ H, Wgt ] = multi_station_observation_weight(currentObservationTime, currentObservationLinkEnds, currentInformationMatrix,...
            variance, 0.5);
        
        % Perform formal error analysis
        inverseNormalizedCovarianceMatrix = H'*Wgt*H + normalizedInverseAprioriCovarianceMatrix;
        
        % Inverse
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
        
        % Take the standard deviation
        sigma0=sqrt(diag(unnormalizedCovarianceMatrix));
        sigmas_comb(:, indexIter, matrixArrayIndex)  = sigma0(7:end);
        indexIter = indexIter + 1;
    end
    
    matrixArrayIndex = matrixArrayIndex + 1;
end

%% Plots

load('formal3.mat')

% The aim is to understand whether the std decrases with the observation time and by including the radio telescopes
figure( 1 )

subplot( 2,1,1 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 1, : ) , ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 1, :, 1 ) , ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 1, :, 2 ) , ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 1, :, 3 ) , ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 1, :, 4 ) , ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{F}$','Interpreter','latex','FontSize',40)
box on
grid on

subplot( 2,1,2 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 2, : ) , ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 2, :, 1 ), ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 2, :, 2 ), ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 2, :, 3 ), ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 2, :, 4 ), ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{FCN}$ [rad/s]','Interpreter','latex','FontSize',40)
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on
set(gcf,'Position', 0.85*get(0, 'Screensize'));
saveas(gcf,'trackstrg1.png')

figure( 2 )

subplot( 3,1,1 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 3, : ) , ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 3, :, 1 ) , ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 3, :, 2 ) , ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 3, :, 3 ) , ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 3, :, 4 ) , ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{x}$ [m]','Interpreter','latex','FontSize',40)
set(gca, 'yscale', 'log')
box on
grid on

subplot( 3,1,2 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 4, : ) , ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 4, :, 1 ) , ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 4, :, 2 ) , ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 4, :, 3 ) , ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 4, :, 4 ) , ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{y}$ [m]','Interpreter','latex','FontSize',40)
set(gca, 'yscale', 'log')
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on

subplot( 3,1,3 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 5, : ) , ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 5, :, 1 ) , ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 5, :, 2 ) , ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 5, :, 3 ) , ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 5, :, 4 ) , ':sb', 'Color', 'c')
xline(365,'--k');
xlabel( 'Observation Time' )
ylabel('$\sigma_{z}$ [m]','Interpreter','latex','FontSize',40)
set(gca, 'yscale', 'log')
box on
grid on

set(gcf,'Position', 0.85*get(0, 'Screensize'));
saveas(gcf,'trackstrg2.png')

figure( 3 )

subplot( 4,2,1 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 6, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 6, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 6, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 6, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 6, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{c1}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

subplot( 4,2,2 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 7, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 7, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 7, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 7, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 7, :, 4 ) / mas, ':sb', 'Color', 'c')
ylabel('$\sigma_{\phi_{s1}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

subplot( 4,2,3 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 8, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 8, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 8, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 8, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 8, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{c2}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

subplot( 4,2,4 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 9, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 9, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 9, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 9, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 9, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{s2}}$ [mas]','Interpreter','latex','FontSize',35)
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on

subplot( 4,2,5 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 10, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 10, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 10, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 10, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 10, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{c3}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

subplot( 4,2,6 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 11, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 11, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 11, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 11, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 11, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{s3}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

subplot( 4,2,7 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 12, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 12, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 12, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 12, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 12, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{c4}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

subplot( 4,2,8 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 13, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 13, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 13, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 13, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 13, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{s4}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

set(gcf,'Position', 0.85*get(0, 'Screensize'));
saveas(gcf,'trackstrg3.png')

figure(4)

subplot( 2,2,1 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 14, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 14, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 14, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 14, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 14, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{XP_{c1}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

subplot( 2,2,2 )
set(gca,'FontSize',20)
hold on 
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 15, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 15, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 15, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 15, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 15, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{XP_{s1}}$ [mas]','Interpreter','latex','FontSize',40)
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on

subplot( 2,2,3 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 16, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 16, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 16, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 16, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 16, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{YP_{c1}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

subplot( 2,2,4 )
set(gca,'FontSize',20)
hold on 
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 17, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 17, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 17, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 17, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 17, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{YP_{s1}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

set(gcf,'Position', 0.85*get(0, 'Screensize'));
saveas(gcf,'trackstrg4.png')

figure(5)

subplot( 2,2,1 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 18, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 18, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 18, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 18, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 18, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{XP_{c2}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

subplot( 2,2,2 )
set(gca,'FontSize',20)
hold on 
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 19, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 19, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 19, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 19, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 19, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{XP_{s2}}$ [mas]','Interpreter','latex','FontSize',40)
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on

subplot( 2,2,3 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 20, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 20, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 20, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 20, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 20, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{YP_{c2}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

subplot( 2,2,4 )
set(gca,'FontSize',20)
hold on 
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 21, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 21, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 21, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 21, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 21, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{YP_{s2}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

set(gcf,'Position', 0.85*get(0, 'Screensize'));
saveas(gcf,'trackstrg5.png')

figure(6)

subplot( 2,2,1 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 21, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 21, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 21, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 21, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 21, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{XC_{c}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

subplot( 2,2,2 )
set(gca,'FontSize',20)
hold on 
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 22, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 22, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 22, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 22, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 22, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{XC_{s}}$ [mas]','Interpreter','latex','FontSize',40)
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on

subplot( 2,2,3 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 23, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 23, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 23, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 23, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 23, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{YC_{c}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

subplot( 2,2,4 )
set(gca,'FontSize',20)
hold on 
plot( evaluationTimes_DSS63 / one_day, sigmas_DSS63( 24, : ) / mas, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 24, :, 1 ) / mas, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 24, :, 2 ) / mas, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 24, :, 3 ) / mas, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / one_day, sigmas_comb( 24, :, 4 ) / mas, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{YC_{s}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

set(gcf,'Position', 0.85*get(0, 'Screensize'));
saveas(gcf,'trackstrg6.png')

% Plot the formal error / apriori constrain against each element
figure(7)
set(gca,'FontSize',20)
hold on
plot( (1:1:length(sigmaAPriori(7:end))) , transpose(sigmas_DSS63( :, end )) ./ sigmaAPriori(7:end), ':sb', 'Color', 'b')
plot( (1:1:length(sigmaAPriori(7:end))) , transpose(sigmas_comb( :, end, 1 )) ./ sigmaAPriori(7:end), ':sb', 'Color', 'k')
plot( (1:1:length(sigmaAPriori(7:end))) , transpose(sigmas_comb( :, end, 2 )) ./ sigmaAPriori(7:end), ':sb', 'Color', 'r')
plot( (1:1:length(sigmaAPriori(7:end))) , transpose(sigmas_comb( :, end, 3 )) ./ sigmaAPriori(7:end), ':sb', 'Color', 'g')
plot( (1:1:length(sigmaAPriori(7:end))) , transpose(sigmas_comb( :, end, 4 )) ./ sigmaAPriori(7:end), ':sb', 'Color', 'c')
ylabel('formal error / a priori constrain')
ax = gca;
ax.XTick = 1:33;
ax.XLim = [0 33];
ax.XTickLabel = {'F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
    '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
    'XC_{c3}', 'XC_{s3}', 'YC_{c3}', 'YC_{s3}','XC_{c4}', 'XC_{s4}', 'YC_{c4}', 'YC_{s4}',...
    'XC_{c5}', 'XC_{s5}', 'YC_{c5}', 'YC_{s5}','Interpreter','latex',};   
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on
set(gcf,'Position', 0.85*get(0, 'Screensize'));
saveas(gcf,'trackstrg7.png')