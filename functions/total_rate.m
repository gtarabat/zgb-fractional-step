% File              : total_rate.m
% Author            : George Arampatzis <garampat@ethz.ch>
% Date              : 27.08.2021
% Last Modified Date: 27.08.2021
% Last Modified By  : George Arampatzis <garampat@ethz.ch>
function r = total_rate(s, i, j, rcd)
%% Description: 
%   Total rate of every reaction that can take place on site (i,j)
%
%% Input:
%   s : lattice, square matrix with values {-1,0,+1} of size NxN
%   i : integer, x-position on lattice
%   j : integer, y-position on lattice
%   rcd : struct containing the reaction rates 
%
%% Output:
%   r :  total rate, sum of all rates
%

r1 = rates_vec(s, i, j, rcd);
r = sum(r1);