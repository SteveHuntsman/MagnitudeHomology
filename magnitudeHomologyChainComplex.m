function [MC,boundary,L] = ...
    magnitudeHomologyChainComplex(d,simplices,lengths,blur,varargin)

% Inputs:
%     d,          a matrix of distances
%     simplices,  a cell of simplices from magnitudeHomologySimplices
%     lengths,    a cell of lengths from magnitudeHomologySimplices
%     blur,       a flag for computing blurred magnitude homology
% Optional input:
%     L,          a tuple of lengths. If this is not supplied, then the
%                 unique lengths of simplices is used. For unweighted
%                 (di)graphs, it's best to avoid this argument, but for
%                 weighted cases it's typically advisable to use this in
%                 conjunction with blurring. This is not handled gracefully
%                 in the code yet.
% Outputs:
%     MC,         a cell of chains
%     boundary,   a cell of boundary maps
%     L,          a tuple of lengths (note that this may also be an input)
%
% Copyright (c) 2023, Steve Huntsman. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

assert(numel(simplices)==numel(lengths),...
    'numel(simplices)==numel(lengths)');
K = numel(simplices)-1;

%% Relation to be applied on lengths and grading
% Use numerically well-behaved equality predicate to avoid any
% floating point badness that might have been inherited in D
if blur % blurred magnitude homology
    rel = @(x1,x2) le(x1,x2);   % x1 <= x2
else    % unblurred magnitude homology
    rel = @(x1,x2) max(x1,x2)-min(x1,x2)<sqrt(eps); % loose ver of x1 == x2    
end

%% Unique lengths
if numel(varargin) == 0
    L = unique(cell2mat(lengths(:))');
elseif numel(varargin) == 1
    L = varargin{1};
else
    warning('superflous inputs: expecting only L as optional; ignoring');
end

%% Chains
MC = cell(K+1,numel(L));
for ell = 1:numel(L)
    if L(ell) == 0
        MC{0+1,ell} = simplices{0+1};
    end
    for k = 0:K
        MC{k+1,ell} = simplices{k+1}(rel(lengths{k+1},L(ell)),:);
    end
end

%% Boundary matrices
boundary = cell(K,numel(L));
for k = 1:K
    for ell = 1:numel(L)
        %% Hack to handle the boundary map for 0-simplices
        if L(ell) == 0 && numel(MC{k,ell}) == 1
            assert(k==1,"k==1");
            boundary{k-1+1,ell} = 0;
            % The second part of this hack just below is performed because
            % [~,b] = homologyReal({0;[];[];[]}) gives b = [0,0,0] while
            % [~,b] = homologyReal({0;0;[];[]}) gives b = [1,0,0]
            if K < 2, warning('K < 2'); end
            boundary{k-1+2,ell} = 0;
            continue;
        end
        %% Initialization
        boundaryInit = sparse(size(MC{k-1+1,ell},1),size(MC{k+1,ell},1));
        boundary{k+1,ell} = boundaryInit;
        %% The actual computation
        for j = 1:(k-1)
            boundary_j = boundaryInit;
            for i = 1:size(MC{k+1,ell},1)
                neg = MC{k+1,ell}(i,j+1-1);
                nil = MC{k+1,ell}(i,j+1);
                pos = MC{k+1,ell}(i,j+1+1);
%                 if rel(d(neg,pos),d(neg,nil)+d(nil,pos))
                if abs(d(neg,pos)-d(neg,nil)-d(nil,pos))<sqrt(eps)
                    del = MC{k+1,ell}(i,setdiff(0:k,j)+1);
                    ind = ismember(MC{k-1+1,ell},del,'rows','legacy');
                    boundary_j(ind,i) = 1; %#ok<SPRIX> 
                end
            end
            boundary{k+1,ell} = boundary{k+1,ell}+((-1)^j)*boundary_j;
        end
        %% Confirm boundary^2 = 0
        if ~isempty(boundary{k,ell}) && ~isempty(boundary{k+1,ell})
            assert(all(boundary{k,ell}*boundary{k+1,ell}==0,'all'),...
                'boundary^2 is nonzero');
        end
    end
end