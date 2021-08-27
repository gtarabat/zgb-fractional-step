% File              : rates_vec.m
% Author            : George Arampatzis <garampat@ethz.ch>
% Date              : 27.08.2021
% Last Modified Date: 27.08.2021
% Last Modified By  : George Arampatzis <garampat@ethz.ch>
function r = rates_vec(s, i, j, rcd)
%% Description: 
%   Rate of every reaction that can take place on site (i,j)
%
%% Input:
%   s : lattice, square matrix with values {-1,0,+1} of size NxN
%   i : integer, x-position on lattice
%   j : integer, y-position on lattice
%   rcd : struct containing the reaction rates 
%
%% Output:
%   r : vector that contains the rates
%

N = size(s,1);

t = [0, 0.25, 0.5, 0.75, 1];

[il, ir, ju, jd] = neighbors_of(i, j, N);

suma1 = (s(i,ju)==0) + (s(i,jd)==0) + (s(il,j)==0) + (s(ir,j)==0);
suma2 = (s(i,ju)<0)  + (s(i,jd)<0)  + (s(il,j)<0)  + (s(ir,j)<0);
suma3 = (s(i,ju)>0)  + (s(i,jd)>0)  + (s(il,j)>0)  + (s(ir,j)>0);

c1 = 1 - abs(s(i,j));
c2 = 0.5 * s(i,j) * rcd.k2;

% CO adsorption
r(1) = c1 *  rcd.k1;

% O2 adsorption
r(2) = c1 * ( 1-rcd.k1 ) * t(suma1+1);

% O + CO -> CO2 and desorption
r(3) = c2 * ( 1+s(i,j) ) * t(suma2+1);

% CO + O -> CO2 and desorption
r(4) = c2 * ( s(i,j)-1 ) * t(suma3+1); 
