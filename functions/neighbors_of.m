% File              : neighbors_of.m
% Author            : George Arampatzis <garampat@ethz.ch>
% Date              : 27.08.2021
% Last Modified Date: 27.08.2021
% Last Modified By  : George Arampatzis <garampat@ethz.ch>
function [il, ir, ju, jd] = neighbors_of(i, j, N)
%% Description: 
%   Neighbors of the (i,j) lattice site with periodic boundary conditions.
%
%% Input:
%   i : integer, x-position on lattice
%   j : integer, y-position on lattice
%   N : integer size of a square lattice NxN
%
%% Output:
%   il : integer, x-position of left neighbor
%   ir : integer, x-position of right neighbor
%   ju : integer, y-position of up neighbor
%   jd : integer, y-position of down neighbor

il = mod(i-1-(N+1),N)+1;
ir = mod(i+1-(N+1),N)+1;

ju = mod(j+1-(N+1),N)+1;
jd = mod(j-1-(N+1),N)+1;
