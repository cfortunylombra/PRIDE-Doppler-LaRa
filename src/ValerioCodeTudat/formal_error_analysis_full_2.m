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
% the number of stations and the observation correlation coefficient. The
% variance-covariance matrix associated to the estimation is plotted too, as
% well as the ratio between the formal error and the a priori
% uncertainties.
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
dataDirectory = strrep(matlab.desktop.editor.getActiveFilename,'src\ValerioCodeTudat\formal_error_analysis_full_2.m','output\POD_LaRa\');

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

% Add the name of the transmitter to the ground station names
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

%% PLOT FORMAL ERRORS

% Create observation times, link-ends and estimation information matrix
obsTimes = [];
obsLinkEnds = [];
estInformationMatrix = [];

% Iterate for each ground station
for groundStationIDindex = 1 : length( groundStationID )
    fprintf('groundStationIDindex:')
    disp(groundStationIDindex)
    
    % Append the obsTimes, obslinkEnds and estInformationMatrix for each ground station
    obsTimes = [obsTimes; observationTimes( observationLinkEnds == groundStationIDindex )];
    obsLinkEnds = [obsLinkEnds; observationLinkEnds( observationLinkEnds == groundStationIDindex )];
    estInformationMatrix = [estInformationMatrix; estimationInformationMatrix( observationLinkEnds == groundStationIDindex, : )];
    
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
        
        % Take the standard deviation
        sigma0=sqrt(diag(unnormalizedCovarianceMatrix));
        
        sigmas_mars_state_comb(:, correlationCoefficientIndex, groundStationIDindex)  = sigma0(1 : 6);

        sigmas_comb(:, correlationCoefficientIndex, groundStationIDindex)  = sigma0(7 : end);  
        
        if groundStationIDindex == length( groundStationID )
            correlationsMatrixArray(:,:,correlationCoefficientIndex) = unnormalizedCovarianceMatrix ./ (sigma0 * sigma0');
        end
    end
end

%% Plot 

load('formal2.mat')

% The variance for each element against the correlation coefficient (performed for each ground station)
col = hsv(length( groundStationID ));
for groundStationIDindex = 1 : length( groundStationID )
    
    figure( 1 ) 
    
    subplot( 2,1,1 ) 
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 1, :, groundStationIDindex ), ':s', 'Color', col( groundStationIDindex, : ))   
    ylabel('$\sigma_{F}$','Interpreter','latex','FontSize',40)
    ylim([0 inf])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,1,2 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 2, :, groundStationIDindex ), ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{FCN}$ [rad/s]','Interpreter','latex','FontSize',40)
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
    '9-station','10-station','11-station','12-station','13-station','14-station','15-station','16-station','Location','best')
    ylim([0 inf])
    xlim([-0.05 0.95])
    box on
    grid on
    set(gcf,'Position', 0.85*get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff1.png')
    
    figure( 2 )
       
    subplot( 3,1,1 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 3, :, groundStationIDindex ), ':s', 'Color', col( groundStationIDindex, : ))   
    ylabel('$\sigma_{x}$ [m]','Interpreter','latex','FontSize',40)
    ylim([0 inf])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 3,1,2 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 4, :, groundStationIDindex ), ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{y}$ [m]','Interpreter','latex','FontSize',40)
    ylim([0 inf])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 3,1,3 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 5, :, groundStationIDindex ), ':s', 'Color', col( groundStationIDindex, : ))
    xlabel( 'Covariance Coefficient' )
    ylabel('$\sigma_{z}$ [m]','Interpreter','latex','FontSize',40)
    ylim([0 inf])
    xlim([-0.05 0.95])
    box on
    grid on
    
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
    '9-station','10-station','11-station','12-station','13-station','14-station','15-station','16-station',...
    'Location',[0.91671875 0.434175200534759 0.0640624999999999 0.168323863636363])
    set(gcf,'Position', 0.85*get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff2.png')
    
    figure( 3 ) 
    
    subplot( 4,2,1 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 6, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{c1}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,2 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 7, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{s1}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,3 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 8, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{c2}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,4 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 9, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{s2}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,5 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 10, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{c3}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,6 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 11, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{s3}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,7 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 12, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{c4}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,8 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 13, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{s4}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
    '9-station','10-station','11-station','12-station','13-station','14-station','15-station','16-station',...
    'Location',[0.0210937411317396 0.444590723735932 0.0630625 0.168323863636363])
    set(gcf,'Position', 0.85*get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff3.png')
    
    figure( 4 )
    
    subplot( 2,2,1 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 14, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XP_{c1}}$ [mas]','Interpreter','latex','FontSize',40)
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,2 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 15, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XP_{s1}}$ [mas]','Interpreter','latex','FontSize',40)
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
    '9-station','10-station','11-station','12-station','13-station','14-station','15-station','16-station','Location','best')
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,3 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 16, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YP_{c1}}$ [mas]','Interpreter','latex','FontSize',40)
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,4 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 17, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YP_{s1}}$ [mas]','Interpreter','latex','FontSize',40)
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    set(gcf,'Position', 0.85*get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff4.png')
        
    figure( 5 )
    
    subplot( 2,2,1 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 18, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XP_{c2}}$ [mas]','Interpreter','latex','FontSize',40)
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,2 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 19, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XP_{s2}}$ [mas]','Interpreter','latex','FontSize',40)
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
    '9-station','10-station','11-station','12-station','13-station','14-station','15-station','16-station','Location','best')
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,3 )
    set(gca,'FontSize',20)
    hold on
    grid on
    
    plot( correlationCoefficient, sigmas_comb( 20, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YP_{c2}}$ [mas]','Interpreter','latex','FontSize',40)
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,4 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 21, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YP_{s2}}$ [mas]','Interpreter','latex','FontSize',40)
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    set(gcf,'Position', 0.85*get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff5.png')
        
    figure( 6 )
    
    subplot( 2,2,1 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 22, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XC_{c3}}$ [mas]','Interpreter','latex','FontSize',40)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,2 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 23, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XC_{s3}}$ [mas]','Interpreter','latex','FontSize',40)
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
    '9-station','10-station','11-station','12-station','13-station','14-station','15-station','16-station','Location','best')
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,3 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 24, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YC_{c3}}$ [mas]','Interpreter','latex','FontSize',40)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,4 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 25, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YC_{s3}}$ [mas]','Interpreter','latex','FontSize',40)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    set(gcf,'Position', 0.85*get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff6.png')
    
    figure( 7 )
    
    subplot( 2,2,1 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 26, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XC_{c4}}$ [mas]','Interpreter','latex','FontSize',40)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,2 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 27, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XC_{s4}}$ [mas]','Interpreter','latex','FontSize',40)
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
    '9-station','10-station','11-station','12-station','13-station','14-station','15-station','16-station','Location','best')
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,3 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 28, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YC_{c4}}$ [mas]','Interpreter','latex','FontSize',40)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,4 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 29, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YC_{s4}}$ [mas]','Interpreter','latex','FontSize',40)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    set(gcf,'Position', 0.85*get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff7.png')
    
    figure( 8 )
    
    subplot( 2,2,1 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 30, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XC_{c5}}$ [mas]','Interpreter','latex','FontSize',40)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,2 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 31, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XC_{s5}}$ [mas]','Interpreter','latex','FontSize',40)
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
    '9-station','10-station','11-station','12-station','13-station','14-station','15-station','16-station','Location','best')
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,3 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 32, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YC_{c5}}$ [mas]','Interpreter','latex','FontSize',40)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,4 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 33, :, groundStationIDindex ) / mas, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YC_{s5}}$ [mas]','Interpreter','latex','FontSize',40)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    set(gcf,'Position', 0.85*get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff8.png')     
end

% Plot formal error / apriori constrain against each element
figure( 9 )
set(gca,'FontSize',20)
hold on
legend_label = ["\rho = 0","\rho = 0.1","\rho = 0.2","\rho = 0.3","\rho = 0.4","\rho = 0.5","\rho = 0.6","\rho = 0.7",...
        "\rho = 0.8","\rho = 0.9"];
for ind = 1 : length( correlationCoefficient )
    plot(1:1:length(sigmaAPriori(7:end)), transpose(sigmas_comb( :, ind, end )) ./ sigmaAPriori(7:end), ':s', 'Color', col( ind, : ))
end
box on
grid on
ylabel('formal error / a priori constrain') 
ax = gca;
ax.XTick = 1:33;
ax.XLim = [0 33];
ax.XTickLabel = {'F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
    '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
    'XC_{c3}', 'XC_{s3}', 'YC_{c3}', 'YC_{s3}','XC_{c4}', 'XC_{s4}', 'YC_{c4}', 'YC_{s4}',...
    'XC_{c5}', 'XC_{s5}', 'YC_{c5}', 'YC_{s5}','Interpreter','latex',};   
legend(legend_label)
set(gcf,'Position', 0.85*get(0, 'Screensize'));
saveas(gcf,'sigma_corrcoeff10.png')


% Plot image of matrix data for correlation = 0.0
figure ( 10 )
set(gca,'FontSize',30)
imagesc( abs( correlationsMatrixArray( :, :, 1) ) );
colorbar
ax = gca;
ax.XTick = 1:39;
ax.XLim = [1 39];
ax.XTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c3}', 'XC_{s3}', 'YC_{c3}', 'YC_{s3}','XC_{c4}', 'XC_{s4}', 'YC_{c4}', 'YC_{s4}',...
        'XC_{c5}', 'XC_{s5}', 'YC_{c5}', 'YC_{s5}', 'Interpreter','latex',};
ax.YTick = 1:39;
ax.YLim = [1 39];
ax.YTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c3}', 'XC_{s3}', 'YC_{c3}', 'YC_{s3}','XC_{c4}', 'XC_{s4}', 'YC_{c4}', 'YC_{s4}',...
        'XC_{c5}', 'XC_{s5}', 'YC_{c5}', 'YC_{s5}', 'Interpreter','latex',};
title('$\rho = 0 $ ','Interpreter','latex')
set(gcf,'Position', 0.85*get(0, 'Screensize'));
saveas(gcf,'cov1_1.png')

% Plot image of matrix data for correlation = 0.5
figure( 11 )
set(gca,'FontSize',20)
imagesc( abs( correlationsMatrixArray( :, :, 6) ) );
colorbar
ax = gca;
ax.XTick = 1:39;
ax.XLim = [1 39];
ax.XTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c3}', 'XC_{s3}', 'YC_{c3}', 'YC_{s3}','XC_{c4}', 'XC_{s4}', 'YC_{c4}', 'YC_{s4}',...
        'XC_{c5}', 'XC_{s5}', 'YC_{c5}', 'YC_{s5}', 'Interpreter','latex',};
ax.YTick = 1:39;
ax.YLim = [1 39];
ax.YTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c3}', 'XC_{s3}', 'YC_{c3}', 'YC_{s3}','XC_{c4}', 'XC_{s4}', 'YC_{c4}', 'YC_{s4}',...
        'XC_{c5}', 'XC_{s5}', 'YC_{c5}', 'YC_{s5}', 'Interpreter','latex',};
title('$\rho = 0.5 $ ','Interpreter','latex')
set(gcf,'Position', 0.85*get(0, 'Screensize'));
saveas(gcf,'cov1_2.png')

% Plot image of matrix data for correlation = 0.9
figure( 12 )
set(gca,'FontSize',20)
imagesc( abs( correlationsMatrixArray( :, :, 10 ) ) );
colorbar
ax = gca;
ax.XTick = 1:39;
ax.XLim = [1 39];
ax.XTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c3}', 'XC_{s3}', 'YC_{c3}', 'YC_{s3}','XC_{c4}', 'XC_{s4}', 'YC_{c4}', 'YC_{s4}',...
        'XC_{c5}', 'XC_{s5}', 'YC_{c5}', 'YC_{s5}', 'Interpreter','latex',};
ax.YTick = 1:39;
ax.YLim = [1 39];
ax.YTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c3}', 'XC_{s3}', 'YC_{c3}', 'YC_{s3}','XC_{c4}', 'XC_{s4}', 'YC_{c4}', 'YC_{s4}',...
        'XC_{c5}', 'XC_{s5}', 'YC_{c5}', 'YC_{s5}', 'Interpreter','latex',};
title('$\rho = 0.9 $ ','Interpreter','latex')
set(gcf,'Position', 0.85*get(0, 'Screensize'));
saveas(gcf,'cov1_3.png')

% Plot image of matrix data - difference correlation between 0.9 and 0.0
figure( 13 )
set(gca,'FontSize',20)
%imagesc((abs(correlationsMatrixArray( :, :, 10 )) - abs(correlationsMatrixArray( :, :, 1 )))); %./abs(correlationsMatrixArray( :, :, 1 )))
imagesc(abs(correlationsMatrixArray( :, :, 10 ) - correlationsMatrixArray( :, :, 1 )));
colorbar 
ax = gca;
ax.XTick = 1:39;
ax.XLim = [1 39];
ax.XTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c3}', 'XC_{s3}', 'YC_{c3}', 'YC_{s3}','XC_{c4}', 'XC_{s4}', 'YC_{c4}', 'YC_{s4}',...
        'XC_{c5}', 'XC_{s5}', 'YC_{c5}', 'YC_{s5}', 'Interpreter','latex',};
ax.YTick = 1:39;
ax.YLim = [1 39];
xtickangle(90)
ax.YTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c3}', 'XC_{s3}', 'YC_{c3}', 'YC_{s3}','XC_{c4}', 'XC_{s4}', 'YC_{c4}', 'YC_{s4}',...
        'XC_{c5}', 'XC_{s5}', 'YC_{c5}', 'YC_{s5}', 'Interpreter','latex',};
set(gcf,'Position', 0.85*get(0, 'Screensize'));
saveas(gcf,'cov2.png')











