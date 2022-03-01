format long; clear; close all; clc;

%% DEFINE DATA DIRECTORY AND LOAD FILES

% Path of all the data
dataDirectory = strrep(matlab.desktop.editor.getActiveFilename,'src\ValerioCodeTudat\formal_error_analysis_test.m','output\POD_LaRa\');

% Load estimation information matrix
estimationInformationMatrix = load( strcat( dataDirectory, 'estimation_information_matrix.dat'));

% Load information matrix transformation diagonal
informationMatrixTransformationDiagonal = load( strcat( dataDirectory,...
    'estimation_information_matrix_normalization.dat'));

% Load overall observation time
observationTimes = load( strcat( dataDirectory, 'concatenated_times.dat' ));
observationTimes = (observationTimes - observationTimes(1) );

% Load observation link ends
observationLinkEnds = load( strcat( dataDirectory, 'concatenated_link_ends.dat' ));
observationLinkEnds = observationLinkEnds + 1;

% List of Ground Station Names
groundStationID = ["DSS-63";"BADARY"; "CEDUNA"; "HARTRAO"; "HART15M"; "HOBART12"; "HOBART26"; "TIANMA65"; "WARK30M"; "EFLSBERG"; "IRBENE"; "YEBES40M"; "MEDICINA"; "WETTZELL"; "ONSALA60"; "WRT0"];

%% Define Inverse Apriopri Covariance 

% 1 day in seconds
one_day = 86400; %seconds

% Define 1 mas
mas = pi / ( 180.0 * 1000.0 * 3600.0 );

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
sigma_1 = 0.05E-3 / c ;
sigma_2 = 0.05E-3 / c ;

% Define possible correlation coefficients
covarianceCoefficient = (0:0.1:0.9);

% Define variance
variance = [sigma_1 sigma_2];

x_Sigmas = {};
y_Sigmas = {};
z_Sigmas = {};
%% Sort Data
for index1 = 1 : length(groundStationID)
    fprintf('groundStationID:')
    disp(index1)
    
    % Take the observation times, estimation information matrix and
    % observation link ends from transmitter
    observationTimes = repmat(observationTimes( observationLinkEnds == 1 ), index1, 1);
    estimationInformationMatrix = repmat(estimationInformationMatrix( observationLinkEnds == 1, : ), index1, 1);
    observationLinkEnds = repmat(observationLinkEnds( observationLinkEnds == 1 ), index1, 1);
    
    [~,idxu,idxc] = unique( observationTimes );
    [count, ~, idxcount] = histcounts( idxc,numel( idxu ) );
    idxkeep = count( idxcount ) == 1;
    
    % Not repeated observation times
    uncommonObservationTimes = observationTimes( idxkeep );
    
    % Create and sort repeated observation times array
    commonObservationTimes = observationTimes;
    commonObservationTimes(idxkeep) = []; % Remove the non-repeated observations
    [commonObservationTimes,sortIndex] = sort( commonObservationTimes );
    
    % Not repeated observation link ends 
    uncommonObservationLinkEnds = observationLinkEnds( idxkeep );
    
    % Repeated observation link ends
    commonObservationLinkEnds = observationLinkEnds( ~idxkeep );
    
    % Not repeated estimation information matrix
    uncommonEstimationInformationMatrix = estimationInformationMatrix( idxkeep, : );
    
    % Sort repeated estimation information matrix
    commonEstimationInformationMatrix = estimationInformationMatrix( ~idxkeep, : );
    commonEstimationInformationMatrix = commonEstimationInformationMatrix( sortIndex, : );
    
    H = [ commonEstimationInformationMatrix ; uncommonEstimationInformationMatrix ];
    
    
    x_sigmas = [];
    y_sigmas = [];
    z_sigmas = [];
    for covarianceCoefficientIndex = 1 : length(covarianceCoefficient)
        fprintf('covarianceCoefficient:')
        disp(covarianceCoefficientIndex)
    
        %% CREATE WEIGHTS MATRIX
        
        % Take only once the repeated observation times
        commonTimes = unique( commonObservationTimes );
        
        % Create a sparse matrix for common observation
        commonObservationWgt = sparse([]);
        
        % Iterate for each repeated observation times
        for currentIndex = 1 : length( commonTimes )
            currentTime = commonTimes( currentIndex );
            
            % Compute the number of repeated observation time
            blockDimension = sum( commonObservationTimes == currentTime );
            
            % Find the indexes for the repeated observation times
            [~, currentCommonObservationLinkEndsIndex] = ismember(currentTime, commonObservationTimes);
            
            % Define the variance array
            currentVarianceArray = [];
            for ind = 0 : blockDimension-1
                % Take the link-end for each repeated observation time
                linkEndId = commonObservationLinkEnds(currentCommonObservationLinkEndsIndex + ind);
                % Concatenate the current variance array
                currentVarianceArray = [ currentVarianceArray; variance(linkEndId) ];
            end
            % Since it is a variance, the array has to be power of 2
            varianceCovarianceBlock = currentVarianceArray * currentVarianceArray';
            
            % Create the covariance matrix
            idx = eye( blockDimension );
            idx = idx + covarianceCoefficient(covarianceCoefficientIndex) * (1-idx);
            varianceCovarianceBlock = idx .* varianceCovarianceBlock;
            
            % Append to the Weighted matrix with the inverse of the covariance matrix
            commonObservationWgt = blkdiag( commonObservationWgt, sparse(inv( varianceCovarianceBlock )));
        end
        
        % If there are not repeated observation times
        if isempty( uncommonObservationTimes ) == 1
            Wgt = commonObservationWgt;
        % If there are repeated observation times
        else
            uncommonObservationWgt = zeros([ length( uncommonObservationTimes ) 1 ]);
            for ind = 1 : length( uncommonObservationTimes )
                uncommonObservationWgt( ind ) = 1 / variance(uncommonObservationLinkEnds(ind))^2 ;
            end
            Wgt = blkdiag( commonObservationWgt, sparse( diag( uncommonObservationWgt )));
        end
     
        % Perform formal error analysis
        inverseNormalizedCovarianceMatrix = H'*Wgt*H +normalizedInverseAprioriCovarianceMatrix; 
        
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

% Plot the formal error ratio
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