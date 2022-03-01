function [ H, Wgt] = ...
    multi_station_observation_weight( observationTimes, observationLinkEnds, estimationInformationMatrix, variance, correlationCoefficient)
%MULTI_STATION_OBSERVATION_WEIGHT - function to weight the observations from different stations at same time (H1 line)
%
% Syntax:  [ H, Wgt] = multi_station_observation_weight( observationTimes, observationLinkEnds, estimationInformationMatrix, variance, correlationCoefficient) 
%
% Outputs:
%    H - Sorted design matrix
%    Wgt - Weight matrix
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Valerio Filice
% email address: filicevalerio@gmail.com  
% Last revision: 20-Oct-2019
%
% Copyright (c) 2019, Valerio Filice - Dominic Dirkx
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

%% ------------- BEGIN CODE --------------
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

% CREATE WEIGHTS MATRIX

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
    idx = idx + correlationCoefficient * (1-idx);
    varianceCovarianceBlock = idx .* varianceCovarianceBlock;
    
    % Append to the Weighted matrix with the inverse of the covariance matrix
    commonObservationWgt = blkdiag( commonObservationWgt, sparse(inv( varianceCovarianceBlock )));
end

% If there are not repeated observation times
if isempty( uncommonObservationTimes ) == 1
    % Weight matrix is:
    Wgt = commonObservationWgt;
% If there are repeated observation times
else
    % Since there are not other station with the same observation times,
    % only the diagonal is added
    uncommonObservationWgt = zeros([ length( uncommonObservationTimes ) 1 ]);
    for ind = 1 : length( uncommonObservationTimes )
        % Diagonal
        uncommonObservationWgt( ind ) = 1 / variance(uncommonObservationLinkEnds(ind))^2 ;
    end
    % Concatenate with the repeated Weighted matrix
    Wgt = blkdiag( commonObservationWgt, sparse( diag( uncommonObservationWgt )));
end