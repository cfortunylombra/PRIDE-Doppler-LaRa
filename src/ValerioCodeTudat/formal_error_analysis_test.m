format long; clear; close all; clc;

%% DEFINE DATA DIRECTORY AND LOAD FILES

dataDirectory = '/Users/valeriofilice/Tudat/tudatBundle/tudatApplications/Outputs/LaRa_orbit_determination11/';

estimationInformationMatrix = load( strcat( dataDirectory, 'EstimationInformationMatrix.dat' ));
informationMatrixTransformationDiagonal = load( strcat( dataDirectory,...
    'EstimationInformationMatrixNormalization.dat' ));
observations = load( strcat( dataDirectory, 'ObservationMeasurements.dat' ));
observationTimes = load( strcat( dataDirectory, 'ObservationTimes.dat' ));
observationTimes = (observationTimes - observationTimes(1) );
formalError = load( strcat( dataDirectory, 'ObservationFormalEstimationError.dat' ));
estimationResiduals = load( strcat( dataDirectory, 'EstimationResiduals.dat' ));
observationLinkEnds = load( strcat( dataDirectory, 'ObservationLinkEnds.dat' ));
observationLinkEnds = observationLinkEnds + 1;
% groundStationID = {'DSS-63', 'TIANMA65', 'EFLSBERG', 'ARECIBO', 'NOTO', 'YEBES', 'ONSALA60',...
%     'FD-VLBA', 'OV-VLBA', 'HOBART26', 'CEDUNA'};
groundStationID = {'DSS-63', 'EFLSBERG'};

%% DEFINE VARIANCE-COVARIANCE MATRIX PARAMETERS

c = 299792458.0;
sigma_1 = 0.05E-3 / c ;
sigma_2 = 0.05E-3 / c ;

% covarianceCoefficient = [0 0.999];
covarianceCoefficient = (0:0.1:0.9);
% covarianceCoefficient = (-0.9:0.1:0.9);
% covarianceCoefficient = (-0.9:0.01:0.9);

variance = [sigma_1 sigma_2];

x_Sigmas = {};
y_Sigmas = {};
z_Sigmas = {};
%% Sort Data
for index1 = 1 : length(groundStationID)
    observationTimes = repmat(observationTimes( observationLinkEnds == 1 ), index1, 1);
    
    estimationInformationMatrix = repmat(estimationInformationMatrix( observationLinkEnds == 1, : ), index1, 1);
    
    observationLinkEnds = repmat(observationLinkEnds( observationLinkEnds == 1 ), index1, 1);
    
    
    [~,idxu,idxc] = unique( observationTimes );
    [count, ~, idxcount] = histcounts( idxc,numel( idxu ) );
    idxkeep = count( idxcount ) == 1;
    
    uncommonObservationTimes = observationTimes( idxkeep );
    
    commonObservationTimes = observationTimes;
    commonObservationTimes(idxkeep) = [];
    [commonObservationTimes,sortIndex] = sort( commonObservationTimes );
    
    uncommonObservationLinkEnds = observationLinkEnds( idxkeep );
    
    commonObservationLinkEnds = observationLinkEnds( ~idxkeep );
    
    uncommonEstimationInformationMatrix = estimationInformationMatrix( idxkeep, : );
    
    commonEstimationInformationMatrix = estimationInformationMatrix( ~idxkeep, : );
    commonEstimationInformationMatrix = commonEstimationInformationMatrix( sortIndex, : );
    
    H = [ commonEstimationInformationMatrix ; uncommonEstimationInformationMatrix ];
    
    
    x_sigmas = [];
    y_sigmas = [];
    z_sigmas = [];
    for covarianceCoefficientIndex = 1 : length(covarianceCoefficient)
        %% CREATE WEIGHTS MATRIX
        
        commonTimes = unique( commonObservationTimes );
        commonObservationWgt = sparse([]);
        for currentIndex = 1 : length( commonTimes )
            
            currentTime = commonTimes( currentIndex );
            
            blockDimension = sum( commonObservationTimes == currentTime );
            
            [~, currentCommonObservationLinkEndsIndex] = ismember(currentTime, commonObservationTimes);
            
            currentVarianceArray = [];
            for ind = 0 : blockDimension-1
                
                linkEndId = commonObservationLinkEnds(currentCommonObservationLinkEndsIndex + ind);
                
                currentVarianceArray = [ currentVarianceArray; variance(linkEndId) ];
                
            end
            
            varianceCovarianceBlock = currentVarianceArray * currentVarianceArray';
            idx = eye( blockDimension );
            idx = idx + covarianceCoefficient(covarianceCoefficientIndex) * (1-idx);
            varianceCovarianceBlock = idx .* varianceCovarianceBlock;
            
            commonObservationWgt = blkdiag( commonObservationWgt, sparse(inv( varianceCovarianceBlock )));
            
        end
        
        if isempty( uncommonObservationTimes ) == 1
            Wgt = commonObservationWgt;
        else
            uncommonObservationWgt = zeros([ length( uncommonObservationTimes ) 1 ]);
            for ind = 1 : length( uncommonObservationTimes )
                uncommonObservationWgt( ind ) = 1 / variance(uncommonObservationLinkEnds(ind))^2 ;
            end
            Wgt = blkdiag( commonObservationWgt, sparse( diag( uncommonObservationWgt )));
        end
        
        %% PERFORM FORMAL ERROR ANALYSIS
        
        inverseNormalizedCovarianceMatrix = H'*Wgt*H ;
        
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
        x_sigma=sigma0(7);
        y_sigma=sigma0(8);
        z_sigma=sigma0(9);
        
        x_sigmas(covarianceCoefficientIndex) = x_sigma;
        y_sigmas(covarianceCoefficientIndex) = y_sigma;
        z_sigmas(covarianceCoefficientIndex) = z_sigma;
                
    end

    x_Sigmas{index1} = x_sigmas;
    y_Sigmas{index1} = y_sigmas;
    z_Sigmas{index1} = z_sigmas;
    
end

fprintf('Formal error ratio \n')
fprintf('sigma_x = \n')
disp((x_Sigmas{2}./x_Sigmas{1}))
fprintf('sigma_y = \n')
disp((y_Sigmas{2}./y_Sigmas{1}))
fprintf('sigma_z = \n')
disp((z_Sigmas{2}./z_Sigmas{1}))

figure(1)
set(gca,'FontSize',20)
hold on
plot(covarianceCoefficient, x_Sigmas{2}./x_Sigmas{1},'o', 'MarkerSize',17)
plot(covarianceCoefficient, y_Sigmas{2}./y_Sigmas{1},'*', 'MarkerSize',17)
plot(covarianceCoefficient, z_Sigmas{2}./z_Sigmas{1},'+', 'MarkerSize',17)
hold off
xlim([-0.05 0.95])
xlabel('Covariance Coefficient')
ylabel('$\frac{\sigma_{two \ stations}}{\sigma_{one \ station}}$','Interpreter','latex','FontSize',40)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
legend('$\sigma_{x}$ ratio', '$\sigma_{y}$ ratio', '$\sigma_{z}$ ratio','Interpreter','latex')
box on
grid on

% '$\sigma_{YC_{c}}$ [mas]','Interpreter','latex','FontSize',40)

