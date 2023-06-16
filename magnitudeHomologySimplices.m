function [simplices,lengths] = magnitudeHomologySimplices(...
    d,Ad,distanceDigraph,s,t,K,maxLength,eulerian)
    
% Inputs:
%     d,          a matrix of distances
%     Ad,         an unweighted adjacency matrix
%     s,          the source vertex
%     t,          the target vertex
%     K,          the maximum dimension for homology
%     maxLength,  the maximum length of simplices
%     eulerian,   a flag for the Eulerian case 
% Outputs:
%     simplices,  a cell of simplices
%     lengths,    a cell of corresponding lengths
%
% Copyright (c) 2023, Steve Huntsman. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

simplices = cell(1,K+1);
lengths = cell(1,K+1);
% Compute number of paths from s to t that have hop lengths from 0 to K,
% i.e., the number of simplices
pow = eye(size(Ad));
num = zeros(1,K+1);
for k = 0:K
    num(k+1) = pow(s,t);
    if k < K
        pow = pow*Ad;
    end
end
%% Get paths of the desired flavor
% Initializing to {[]} lets us reuse code for simplices and lengths below.
paths = {[]};
if ~eulerian
    if sum(num) > 0
        % NB. Using allpaths and allcycles would be tricky and burdensome
        % here, though the former is actually exceptionally convenient for
        % the Eulerian case just below. Fortunately we can freely use
        % reashortestpaths from
        % https://github.com/SteveHuntsmanBAESystems/BasicPathHomology/.
        % Since reashortestpaths returns a row cell, we'll follow with
        %   paths = paths(:); 
        paths = reashortestpaths(distanceDigraph,s,t,sum(num));
    end
else
    paths = allpaths(distanceDigraph,s,t,'MaxPathLength',maxLength);
end
paths = paths(:);   % to enable cell2mat 

%% Simplices and lengths
hops = cellfun(@numel,paths)-1;
for k = 0:K
    simplices{k+1} = cell2mat(paths(hops==k));
    lengths{k+1} = zeros(size(simplices{k+1},1),1);
    for j = 1:size(simplices{k+1},1)
        for i = 1:(size(simplices{k+1},2)-1)
            lengths{k+1}(j) = lengths{k+1}(j)...
                +d(simplices{k+1}(j,i),simplices{k+1}(j,i+1));
        end
    end
end