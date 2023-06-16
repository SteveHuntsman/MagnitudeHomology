function [allReps,betti,L,allSimplices,allLengths,allMC,allBoundary] = ....
    magnitudeHomology(D,K,blur,eulerian,varargin)

% Magnitude homology calculations (over the reals) for directed weighted
% graphs. Includes flags to compute blurred and Eulerian magnitude
% homology. Aggregates Betti numbers over all sources and targets, but
% produces representatives for individual sources and targets.
%
% Inputs:
%     D,          a (di)graph, possibly arc-weighted (by length)
%     K,          the maximum dimension for homology
%     blur,       a flag for computing blurred magnitude homology
%     eulerian,   a flag for the Eulerian case 
% Optional input:
%     L,          a tuple of lengths. If this is not supplied, then the
%                 unique lengths of simplices is used. For unweighted
%                 (di)graphs, it's best to avoid this argument, but for
%                 weighted cases it's typically advisable to use this in
%                 conjunction with blurring. 
% Outputs:
%     allReps,    a cell of homology representatives, indexed as 
%                 allReps{s,t}{ell,k+1}, where s, t, ell, k+1 respectively
%                 indicate source, target, length, and dimension (plus 1)
%     betti,      an array of Betti numbers
%     L,          the unique finite lengths, or L passed through if an arg
%     allSimplices,
%                 a cell array of simplices, indexed by source and target,
%                 and with each entry itself a cell of simplices of lengths
%                 specified by L
%     allLengths, similarly for lengths of simplices
%     allMC,      a cell array of simplices, indexed by source and target,
%                 containing chains
%     allBoundary,
%                 similarly for boundary maps
%
% % Example for a DAG:
% directsum = @(A,B) [A,zeros(size(A,1),size(B,2));zeros(size(B,1),size(A,2)),B];
% foo = directsum(directsum(ones(5,4),ones(4,3)),ones(3,2));
% bar = [zeros(size(foo,1),5),foo;zeros(2,size(foo,2)+5)];
% D = digraph(bar);
% [allReps,betti,L,allSimplices,allLengths,allMC,allBoundary] = ...
%     magnitudeHomology(D,5,0,0,0:6); betti
% % yields betti = 
% %     14     0     0     0     0     0     0
% %      0    38     0     0     0     0     0
% %      0     0    61     0     0     0     0
% %      0     0     0    60     0     0     0
% %      0     0     0     0     0     0     0
% % in agreement with a prior calculation of the author using different
% % code: see https://hal.science/hal-03688966
%
% % Example for a strong digraph:
% D = digraph([0,1,1;1,0,1;1,0,0]);
% [allReps,betti,L,allSimplices,allLengths,allMC,allBoundary] = ...
%     magnitudeHomology(D,5,0,0,0:6); betti
% % yields betti = 
% %      3     0     0     0     0     0     0
% %      0     5     0     0     0     0     0
% %      0     0     7     0     0     0     0
% %      0     0     0     9     0     0     0
% %      0     0     0     0    11     0     0
% % in agreement with a calculation of D. Govc: see slide 15 of
% % https://mat1.uibk.ac.at/heiko/mag19/Govc.pdf
%
% % Example for Eulerian magnitude homology:
% D = graph([0,1,2,3]+1,[1,2,3,1]+1);
% [allReps,betti,L,allSimplices,allLengths,allMC,allBoundary] = ...
%     magnitudeHomology(D,5,0,1,0:6); betti
% % yields betti = 
% %      4     0     0     0     0     0     0
% %      0     8     0     0     0     0     0
% %      0     0     6     8     0     0     0
% %      0     0     0     0    14     4     0
% %      0     0     0     0     0     0     0
% % in agreement with a calculation of G. Menara: see
% % https://golem.ph.utexas.edu/category/2023/04/eulerian_magnitude_homology.html
%
% Copyright (c) 2023, Steve Huntsman. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

%% If D has no weights, add unit weights; remove any loops and get size
% Weights
foo = string(D.Edges.Properties.VariableNames);
if nnz(foo=="Weight") == 0
    D.Edges.Weight = ones(size(D.Edges,1),1);
end
% Loops
D = simplify(D);

%% Get distance-weighted reachability relation in digraph & adjacency form
d = distances(D);
% Reachability weighted by distances
distanceDigraph = digraph(d);
distanceDigraph = ...
    rmedge(distanceDigraph,find(isinf(distanceDigraph.Edges.Weight)));
Ad = adjacency(distanceDigraph); % NOT weighted

%% Handle length grading and blurring
if numel(varargin) == 0
    L = unique(d(~isinf(d)));
elseif numel(varargin) == 1
    L = varargin{1};
    if blur == 0
        ud = unique(d); 
        if any(ud~=round(ud))
            warning('distances not all integers: blurring is a good idea');
        end
        if numel(setdiff(ud,L)) > 0
            warning('distances not all in L');
        end
    end
else
    warning('superflous inputs: expecting only L as optional; ignoring');
end
maxLength = max(L);

%% Loop over sources and targets for simplices, chains, boundary maps, etc.
allReps = cell(size(D.Nodes,1),size(D.Nodes,1)); 
betti = zeros(K,numel(L));
allSimplices = cell(size(D.Nodes,1),size(D.Nodes,1));
allLengths = cell(size(D.Nodes,1),size(D.Nodes,1));
allMC = cell(size(D.Nodes,1),size(D.Nodes,1));
allBoundary = cell(size(D.Nodes,1),size(D.Nodes,1));
for s = 1:size(D.Nodes,1)
    srcString = ['src = ',num2str(s),'/',num2str(size(D.Nodes,1))];
    for t = 1:size(D.Nodes,1)
        tarString = ['tar = ',num2str(t),'/',num2str(size(D.Nodes,1))];
        disp([srcString,'; ',tarString]);
        [simplices,lengths] = magnitudeHomologySimplices(...
            d,Ad,distanceDigraph,s,t,K,maxLength,eulerian);
        [MC,boundary,~] = magnitudeHomologyChainComplex(...
            d,simplices,lengths,blur,L); % L is input: don't assign 3rd arg
        %%
        allSimplices{s,t} = simplices;
        allLengths{s,t} = lengths;
        allMC{s,t} = MC;
        allBoundary{s,t} = boundary;
        %% 
        betti_st = zeros(K,numel(L));
        for ell = 1:numel(L)
            [hom_ell,betti_ell] = homologyReal(boundary(:,ell));
            betti_st(:,ell) = betti_ell;
            allReps{s,t}(:,ell) = hom_ell;
        end
        betti = betti+betti_st;
    end
end