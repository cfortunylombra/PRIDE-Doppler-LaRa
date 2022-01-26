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

%% Inverse Apriopri Covariance

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

% inverseAPrioriCovariance = ...
%     diag([1 1 1 1 1 1 ...
%     1/(0.07)^2 1/(sigma_FCN_0)^2 1/(200E3)^2 1/(200E3)^2 1/(200E3)^2 ...
%     1/(milliArcSecondToRadian * 10)^2 1/(milliArcSecondToRadian * 12)^2 1/(milliArcSecondToRadian * 9)^2 ...
%     1/(milliArcSecondToRadian * 8)^2 1/(milliArcSecondToRadian * 8)^2 1/(milliArcSecondToRadian * 7)^2 ...
%     1/(milliArcSecondToRadian * 6)^2 1/(milliArcSecondToRadian * 6)^2 ...
%     1/(2.8 * milliArcSecondToRadian * sin( deg2rad( 46.5 )))^2 1/(2.8 * milliArcSecondToRadian * cos( deg2rad( 46.5 )))^2 ...
%     1/(11.7 * milliArcSecondToRadian * sin( deg2rad( 118.7 )))^2 1/(11.7 * milliArcSecondToRadian * cos( deg2rad( 118.7 )))^2 ...
%     1/(8.9 * milliArcSecondToRadian * sin( deg2rad( -150.1 )))^2 1/(11.7 * milliArcSecondToRadian * sin( deg2rad( -150.1 )))^2 ...
%     1/(3.9 * milliArcSecondToRadian * sin( deg2rad( 172.5 )))^2 1/(11.7 * milliArcSecondToRadian * sin( deg2rad( 172.5 )))^2 ...
%     1/(50*milliArcSecondToRadian)^2 1/(50*milliArcSecondToRadian)^2 ...
%     1/(50*milliArcSecondToRadian)^2 1/(50*milliArcSecondToRadian)^2]);  


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

%% PLOT FORMAL ERRORS - TRANSMITTER OBSERVATIONS ONLY / PRIDE OBSERVATIONS ONLY / TRANSMITTER + PRIDE OBSERVATIONS

step = divisors( length( transmitterObservationTimes ));
step = step( end - 2 );

evaluationTimes_DSS63 = [];

transmitterWeightMatrix = eye( length( transmitterObservationTimes )) * 1 / sigma_DSS63^2;

indexIter = 1;
for Set = 0 : step : length( transmitterObservationTimes )
    
    if Set == length( transmitterObservationTimes )
        
        currentObservationTime = transmitterObservationTimes( end );
        
        currentInformationMatrix = transmitterInformationMatrix; 
    
        currentWeightMatrix = transmitterWeightMatrix;
        
    elseif Set == 0
        
        currentObservationTime = 0;
    
        currentInformationMatrix = 0;
        
        currentWeightMatrix = 0;
               
    else
        
        currentObservationTime = transmitterObservationTimes( Set );
    
        currentInformationMatrix = transmitterInformationMatrix( 1 : (Set + step - 1), : );
        
        currentWeightMatrix = transmitterWeightMatrix( 1 : (Set + step - 1), 1 : (Set + step - 1) );
        
    end
        
    inverseNormalizedCovarianceMatrix = currentInformationMatrix' * currentWeightMatrix * currentInformationMatrix + ...
        normalizedInverseAprioriCovarianceMatrix;
    
    normalizedCovarianceMatrix = inv(inverseNormalizedCovarianceMatrix);
    
    unnormalizedCovarianceMatrix = zeros(size(normalizedCovarianceMatrix));
    for index=1:length(informationMatrixTransformationDiagonal)
        
        for ind=1:length(informationMatrixTransformationDiagonal)
            
            unnormalizedCovarianceMatrix(index,ind) = ...
                normalizedCovarianceMatrix(index,ind) / (informationMatrixTransformationDiagonal(index) * ...
                informationMatrixTransformationDiagonal(ind));
        end
        
    end
    
    sigma0 = sqrt(diag(unnormalizedCovarianceMatrix))';
    
    sigmas_DSS63(:, indexIter)  = sigma0(7:end);
       
    evaluationTimes_DSS63 = [ evaluationTimes_DSS63 currentObservationTime ];

    indexIter = indexIter + 1;
end

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
    
    indexIter = 1;
    
    for Set = 1 : length( evaluationTimes_DSS63 )
        
        if Set == length( evaluationTimes_DSS63 )
            
            currentObservationTime = time;
            
            currentInformationMatrix = informationMatrix;
            
            currentObservationLinkEnds =  linkEnds;
            
        else
            
            currentObservationTime = time( time <= evaluationTimes_DSS63( Set ));
            
            currentInformationMatrix = informationMatrix( time <= evaluationTimes_DSS63( Set ), : );
            
            currentObservationLinkEnds = linkEnds( time <= evaluationTimes_DSS63( Set ));
            
            
        end
        
        [ H, Wgt ] = multi_station_observation_weight( currentObservationTime, currentObservationLinkEnds, currentInformationMatrix,...
            variance, 0.5);
        
        
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
        
        sigma0=sqrt(diag(unnormalizedCovarianceMatrix));

        sigmas_comb(:, indexIter, matrixArrayIndex)  = sigma0(7:end);

        indexIter = indexIter + 1;
        
    end
    
    matrixArrayIndex = matrixArrayIndex + 1;
        
end

%% Plots

figure( 1 )

subplot( 2,1,1 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 1, : ) , ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 1, :, 1 ) , ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 1, :, 2 ) , ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 1, :, 3 ) , ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 1, :, 4 ) , ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{F}$','Interpreter','latex','FontSize',40)
box on

subplot( 2,1,2 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 2, : ) * 2*pi/FCN_0^2/86400, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 2, :, 1 ) * 2*pi/FCN_0^2/86400, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 2, :, 2 ) * 2*pi/FCN_0^2/86400, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 2, :, 3 ) * 2*pi/FCN_0^2/86400, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 2, :, 4 ) * 2*pi/FCN_0^2/86400, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{FCN}$ [days]','Interpreter','latex','FontSize',40)
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'trackstrg1.png')

figure( 2 )

subplot( 3,1,1 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 3, : ) , ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 3, :, 1 ) , ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 3, :, 2 ) , ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 3, :, 3 ) , ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 3, :, 4 ) , ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{x}$ [m]','Interpreter','latex','FontSize',40)
set(gca, 'yscale', 'log')
box on
grid on

subplot( 3,1,2 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 4, : ) , ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 4, :, 1 ) , ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 4, :, 2 ) , ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 4, :, 3 ) , ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 4, :, 4 ) , ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{y}$ [m]','Interpreter','latex','FontSize',40)
set(gca, 'yscale', 'log')
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on

subplot( 3,1,3 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 5, : ) , ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 5, :, 1 ) , ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 5, :, 2 ) , ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 5, :, 3 ) , ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 5, :, 4 ) , ':sb', 'Color', 'c')
xline(365,'--k');
xlabel( 'Observation Time' )
ylabel('$\sigma_{z}$ [m]','Interpreter','latex','FontSize',40)
set(gca, 'yscale', 'log')
box on
grid on

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'trackstrg2.png')

figure( 3 )

subplot( 4,2,1 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 6, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 6, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 6, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 6, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 6, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{c1}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

subplot( 4,2,2 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 7, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 7, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 7, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 7, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 7, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
ylabel('$\sigma_{\phi_{s1}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

subplot( 4,2,3 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 8, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 8, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 8, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 8, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 8, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{c2}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

subplot( 4,2,4 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 9, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 9, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 9, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 9, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 9, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{s2}}$ [mas]','Interpreter','latex','FontSize',35)
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on

subplot( 4,2,5 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 10, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 10, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 10, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 10, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 10, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{c3}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

subplot( 4,2,6 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 11, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 11, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 11, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 11, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 11, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{s3}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

subplot( 4,2,7 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 12, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 12, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 12, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 12, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 12, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{c4}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

subplot( 4,2,8 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 13, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 13, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 13, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 13, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 13, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{\phi_{s4}}$ [mas]','Interpreter','latex','FontSize',35)
box on
grid on

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'trackstrg3.png')

figure(4)

subplot( 2,2,1 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 14, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 14, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 14, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 14, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 14, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{XP_{c1}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

subplot( 2,2,2 )
set(gca,'FontSize',20)
hold on 
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 15, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 15, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 15, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 15, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 15, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{XP_{s1}}$ [mas]','Interpreter','latex','FontSize',40)
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on

subplot( 2,2,3 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 16, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 16, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 16, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 16, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 16, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{YP_{c1}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

subplot( 2,2,4 )
set(gca,'FontSize',20)
hold on 
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 17, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 17, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 17, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 17, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 17, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{YP_{s1}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'trackstrg4.png')

figure(5)

subplot( 2,2,1 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 18, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 18, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 18, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 18, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 18, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{XP_{c2}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

subplot( 2,2,2 )
set(gca,'FontSize',20)
hold on 
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 19, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 19, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 19, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 19, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 19, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{XP_{s2}}$ [mas]','Interpreter','latex','FontSize',40)
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on

subplot( 2,2,3 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 20, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 20, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 20, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 20, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 20, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{YP_{c2}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

subplot( 2,2,4 )
set(gca,'FontSize',20)
hold on 
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 21, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 21, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 21, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 21, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 21, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{YP_{s2}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'trackstrg5.png')

figure(6)

subplot( 2,2,1 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 21, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 21, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 21, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 21, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 21, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{XC_{c}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

subplot( 2,2,2 )
set(gca,'FontSize',20)
hold on 
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 22, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 22, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 22, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 22, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 22, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{XC_{s}}$ [mas]','Interpreter','latex','FontSize',40)
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on

subplot( 2,2,3 )
set(gca,'FontSize',20)
hold on
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 23, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 23, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 23, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 23, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 23, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{YC_{c}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

subplot( 2,2,4 )
set(gca,'FontSize',20)
hold on 
plot( evaluationTimes_DSS63 / 86400, sigmas_DSS63( 24, : ) / milliArcSecondToRadian, ':sb', 'Color', 'b')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 24, :, 1 ) / milliArcSecondToRadian, ':sb', 'Color', 'k')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 24, :, 2 ) / milliArcSecondToRadian, ':sb', 'Color', 'r')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 24, :, 3 ) / milliArcSecondToRadian, ':sb', 'Color', 'g')
plot( evaluationTimes_DSS63 / 86400, sigmas_comb( 24, :, 4 ) / milliArcSecondToRadian, ':sb', 'Color', 'c')
xline(365,'--k');
ylabel('$\sigma_{YC_{s}}$ [mas]','Interpreter','latex','FontSize',40)
box on
grid on

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'trackstrg6.png')

figure(7)
set(gca,'FontSize',20)
hold on
plot( (1:1:length(sigmaAPriori)) , sigmas_DSS63( :, end ) ./ sigmaAPriori, ':sb', 'Color', 'b')
plot( (1:1:length(sigmaAPriori)) , sigmas_comb( :, end, 1 ) ./ sigmaAPriori, ':sb', 'Color', 'k')
plot( (1:1:length(sigmaAPriori)) , sigmas_comb( :, end, 2 ) ./ sigmaAPriori, ':sb', 'Color', 'r')
plot( (1:1:length(sigmaAPriori)) , sigmas_comb( :, end, 3 ) ./ sigmaAPriori, ':sb', 'Color', 'g')
plot( (1:1:length(sigmaAPriori)) , sigmas_comb( :, end, 4 ) ./ sigmaAPriori, ':sb', 'Color', 'c')
ylabel('formal error / a priori constrain')
ax = gca;
ax.XTick = 1:25;
ax.XLim = [0 26];
ax.XTickLabel = {'F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
    '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
    'XC_{c}', 'XC_{s}', 'YC_{c}', 'YC_{s}', 'Interpreter','latex',};
legend('DSS63','DSS63-PRIDE 1w','DSS63-PRIDE 2w','DSS63-PRIDE 3w','DSS63-PRIDE 4w')
box on
grid on
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'trackstrg7.png')