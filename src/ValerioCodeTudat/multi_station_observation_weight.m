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

% CREATE WEIGHTS MATRIX

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
    idx = idx + correlationCoefficient * (1-idx);
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