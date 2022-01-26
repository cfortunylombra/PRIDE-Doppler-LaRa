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

% groundStationID = {'DSS-63', 'EFLSBERG'};

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

%% PLOT FORMAL ERRORS

obsTimes = [];
obsLinkEnds = [];
estInformationMatrix = [];

for groundStationIDindex = 1 : length( groundStationID )

    obsTimes = [obsTimes; observationTimes( observationLinkEnds == groundStationIDindex )];
    
    obsLinkEnds = [obsLinkEnds; observationLinkEnds( observationLinkEnds == groundStationIDindex )];

    estInformationMatrix = [estInformationMatrix; estimationInformationMatrix( observationLinkEnds == groundStationIDindex, : )];
    
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
        
        sigma0=sqrt(diag(unnormalizedCovarianceMatrix));
        
        sigmas_mars_state_comb(:, correlationCoefficientIndex, groundStationIDindex)  = sigma0(1 : 6);

        sigmas_comb(:, correlationCoefficientIndex, groundStationIDindex)  = sigma0(7 : end);  
        
        if groundStationIDindex == length( groundStationID )
            
            correlationsMatrixArray(:,:,correlationCoefficientIndex) = unnormalizedCovarianceMatrix ./ (sigma0 * sigma0');
            
        end
    
    end
            
end

%%

col = hsv(length( groundStationID ));
for groundStationIDindex = 1 : length( groundStationID )
    
    figure( 1 ) 
    
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 1, :, groundStationIDindex ), ':s', 'Color', col( groundStationIDindex, : )) 
    xlabel('Correlation coefficient $\rho_{ij}$','Interpreter','latex')
    ylabel('$\sigma_{F}$','Interpreter','latex','FontSize',40)
    ylim([0 0.045])
    xlim([-0.05 0.95])
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
    '9-station','10-station','11-station','Location','best')
    box on
    grid on
    set(gcf, 'Position', get(0, 'Screensize'))
end

%% Plot

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
    plot( correlationCoefficient, sigmas_comb( 2, :, groundStationIDindex ) * 2*pi/FCN_0^2/86400, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{FCN}$ [days]','Interpreter','latex','FontSize',40)
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
    '9-station','10-station','11-station','Location','best')
    ylim([0 inf])
    xlim([-0.05 0.95])
    box on
    grid on
    set(gcf, 'Position', get(0, 'Screensize'));
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
    '9-station','10-station','11-station','Position',[0.91671875 0.434175200534759 0.0640624999999999 0.168323863636363])
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff2.png')
    
    figure( 3 ) 
    
    subplot( 4,2,1 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 6, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{c1}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,2 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 7, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{s1}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,3 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 8, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{c2}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,4 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 9, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{s2}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,5 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 10, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{c3}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,6 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 11, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{s3}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,7 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 12, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{c4}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 4,2,8 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 13, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{\phi_{s4}}$ [mas]','Interpreter','latex','FontSize',30)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
        '9-station','10-station','11-station','Position',[0.0210937411317396 0.444590723735932 0.0630625 0.168323863636363])
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff3.png')
    
    figure( 4 )
    
    subplot( 2,2,1 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 14, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XP_{c1}}$ [mas]','Interpreter','latex','FontSize',40)
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,2 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 15, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XP_{s1}}$ [mas]','Interpreter','latex','FontSize',40)
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
    '9-station','10-station','11-station','Location','best')
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,3 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 16, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YP_{c1}}$ [mas]','Interpreter','latex','FontSize',40)
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,4 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 17, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YP_{s1}}$ [mas]','Interpreter','latex','FontSize',40)
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff4.png')
        
    figure( 5 )
    
    subplot( 2,2,1 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 18, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XP_{c2}}$ [mas]','Interpreter','latex','FontSize',40)
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,2 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 19, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XP_{s2}}$ [mas]','Interpreter','latex','FontSize',40)
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
    '9-station','10-station','11-station','Location','best')
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,3 )
    set(gca,'FontSize',20)
    hold on
    grid on
    
    plot( correlationCoefficient, sigmas_comb( 20, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YP_{c2}}$ [mas]','Interpreter','latex','FontSize',40)
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,4 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 21, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YP_{s2}}$ [mas]','Interpreter','latex','FontSize',40)
%     ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff5.png')
        
    figure( 6 )
    
    subplot( 2,2,1 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 22, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XC_{c}}$ [mas]','Interpreter','latex','FontSize',40)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,2 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 23, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{XC_{s}}$ [mas]','Interpreter','latex','FontSize',40)
    legend('1-station','2-station','3-station','4-station','5-station','6-station','7-station','8-station',...
    '9-station','10-station','11-station','Location','best')
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,3 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 24, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YC_{c}}$ [mas]','Interpreter','latex','FontSize',40)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    
    subplot( 2,2,4 )
    set(gca,'FontSize',20)
    hold on
    plot( correlationCoefficient, sigmas_comb( 25, :, groundStationIDindex ) / milliArcSecondToRadian, ':s', 'Color', col( groundStationIDindex, : ))
    ylabel('$\sigma_{YC_{s}}$ [mas]','Interpreter','latex','FontSize',40)
    ylim([ 0 inf ])
    xlim([-0.05 0.95])
    box on
    grid on
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff6.png')
            
end

figure( 7 )
set(gca,'FontSize',20)
hold on
for ind = 1 : length( correlationCoefficient )
    
    plot(1:1:length(sigmaAPriori), sigmas_comb( :, ind, end ) ./ sigmaAPriori, ':s', 'Color', col( ind, : ))
    ylabel('formal error / a priori constrain')
    ax = gca;
    ax.XTick = 1:25;
    ax.XLim = [0 26];
    ax.XTickLabel = {'F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c}', 'XC_{s}', 'YC_{c}', 'YC_{s}', 'Interpreter','latex',};
    
    legend('\rho = 0','\rho = 0.1','\rho = 0.2','\rho = 0.3','\rho = 0.4','\rho = 0.5','\rho = 0.6','\rho = 0.7',...
        '\rho = 0.8','\rho = 0.9')
    box on
    grid on
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,'sigma_corrcoeff7.png')
    
end

figure ( 8 )
set(gca,'FontSize',30)
imagesc( abs( correlationsMatrixArray( :, :, 1) ) );
colorbar
ax = gca;
ax.XTick = 1:31;
ax.XLim = [1 31];
ax.XTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c}', 'XC_{s}', 'YC_{c}', 'YC_{s}', 'Interpreter','latex',};
ax.YTick = 1:31;
ax.YLim = [1 31];
ax.YTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c}', 'XC_{s}', 'YC_{c}', 'YC_{s}', 'Interpreter','latex',};
title('$\rho = 0 $ ','Interpreter','latex')
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'cov1_1.png')

figure( 9 )
set(gca,'FontSize',20)
imagesc( abs( correlationsMatrixArray( :, :, 6) ) );
colorbar
ax = gca;
ax.XTick = 1:31;
ax.XLim = [1 31];
ax.XTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c}', 'XC_{s}', 'YC_{c}', 'YC_{s}', 'Interpreter','latex',};
ax.YTick = 1:31;
ax.YLim = [1 31];
ax.YTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c}', 'XC_{s}', 'YC_{c}', 'YC_{s}', 'Interpreter','latex',};
title('$\rho = 0.5 $ ','Interpreter','latex')
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'cov1_2.png')

figure( 10 )
set(gca,'FontSize',20)
imagesc( abs( correlationsMatrixArray( :, :, 10 ) ) );
colorbar
ax = gca;
ax.XTick = 1:31;
ax.XLim = [1 31];
ax.XTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c}', 'XC_{s}', 'YC_{c}', 'YC_{s}', 'Interpreter','latex',};
ax.YTick = 1:31;
ax.YLim = [1 31];
ax.YTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c}', 'XC_{s}', 'YC_{c}', 'YC_{s}', 'Interpreter','latex',};
title('$\rho = 0.9 $ ','Interpreter','latex')
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'cov1_3.png')

figure( 11 )
set(gca,'FontSize',20)
%imagesc((abs(correlationsMatrixArray( :, :, 10 )) - abs(correlationsMatrixArray( :, :, 1 )))); %./abs(correlationsMatrixArray( :, :, 1 )))
imagesc(abs(correlationsMatrixArray( :, :, 10 ) - correlationsMatrixArray( :, :, 1 )));
colorbar 
ax = gca;
ax.XTick = 1:31;
ax.XLim = [1 31];
ax.XTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c}', 'XC_{s}', 'YC_{c}', 'YC_{s}', 'Interpreter','latex',};
ax.YTick = 1:31;
ax.YLim = [1 31];
xtickangle(90)
ax.YTickLabel = {'x_{Mars}','y_{Mars}','z_{Mars}','Vx_{Mars}','Vy_{Mars}','Vz_{Mars}','F','FCN','X_l','Y_l','Z_l','\phi_{c1}', '\phi_{s1}','\phi_{c2}','\phi_{s2}','\phi_{c3}','\phi_{s3}',...
        '\phi_{c4}','\phi_{s4}','XP_{c1}', 'XP_{s1}', 'YP_{c1}', 'YP_{s1}', 'XP_{c2}', 'XP_{s2}', 'YP_{c2}', 'YP_{s2}'...
        'XC_{c}', 'XC_{s}', 'YC_{c}', 'YC_{s}', 'Interpreter','latex',};
saveas(gcf,'cov2.png')
        
set(gcf, 'Position', get(0, 'Screensize'));











