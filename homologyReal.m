function [rep,betti] = homologyReal(boundary)

% Homology over the reals. Uses numerics; can give misleading results. A
% more robust symbolic variant is left as an exercise.
%
% Input:
%     boundary,   a cell of boundary maps
% Outputs:
%     rep,        a cell of homology representatives (over $\mathbb{R}$)
%     betti,      an array of Betti numbers (over $\mathbb{R}$)
%
% Copyright (c) 2023, Steve Huntsman. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

%% 
K = numel(boundary)-1;
ker = cellfun(@(x)null(x,'r'),boundary,'UniformOutput',false);
rep = cell(1,K);
for k = 0:(K-1)
    try
        rep{k+1} = ker{k+1}*null(boundary{k+2}'*ker{k+1},'r');
    catch
        rep{k+1} = [];
    end
end
betti = cellfun(@(x)size(x,2),rep);