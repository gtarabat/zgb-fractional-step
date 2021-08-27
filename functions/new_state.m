% File              : new_state.m
% Author            : George Arampatzis <garampat@ethz.ch>
% Date              : 27.08.2021
% Last Modified Date: 27.08.2021
% Last Modified By  : George Arampatzis <garampat@ethz.ch>
function [ri, rj, s, delta_pop] = new_state(s, i, j, rcd)
%% Description: 
%   Given a site (i,j) of a lattice with size NxN, chooses one of the
%   possible reactions. The function writes directly on the matrix that
%   contains the lattice.
%
%% Input:
%   s : lattice, square matrix with values {-1,0,+1} of size NxN
%   i : integer, x-position on lattice
%   j : integer, y-position on lattice
%   rcd : struct containing the reaction rates 
%
%% Output:
%   ri : x-position of the neighbor site, if the reactions affects a
%        neighbor sitte, -1 otherwise
%   ri : y-position of the neighbor site, if the reactions affects a
%        neighbor sitte, -1 otherwise
%   s  : lattice after the reaction takes place
%   delta_pop : change in the population of every species after the
%               reaction takes place

N = size(s,1);

t = [0, 0.25, 0.5, 0.75, 1];

[il, ir, ju, jd] = neighbors_of(i, j, N);
ii = [ i,  i,  il, ir ]; 
jj = [ jd, ju, j,  j ];

% empty site
if( s(i,j)==0 )
    % find emty neighbor sites
    suma = (s(i,ju)==0) + (s(i,jd)==0) + (s(il,j)==0) + (s(ir,j)==0);
    
    c1 = rcd.k1;
    c2 = (1 - rcd.k1) * t(suma+1);
    c1 = c1 / ( c1 + c2 );
    
    % CO adsorption
    if( rand < c1 )
        s(i,j) = 1;
        ri = -1; 
        rj = -1;
        delta_pop = [ 0, -1, 1 ];
    
    % O2 adsorption
    else
        r = ceil(suma*rand);
        cnt = 0;
        k = 0;
        while( cnt ~= r )
            k = k+1;
            if( s(ii(k),jj(k)) == 0 )
                cnt = cnt+1;
            end
        end
        s(i,j) = -1;
        s( ii(k), jj(k) ) = -1;
        ri = ii(k); 
        rj = jj(k);
        delta_pop = [ 2, -2, 0 ];
    end
    
% CO site: CO + O --> CO2 + desorption
elseif( s(i,j) == 1 )
    suma = (s(i,ju)<0)  + (s(i,jd)<0)  + (s(il,j)<0)  + (s(ir,j)<0);
    
    r = ceil(suma*rand);
    cnt = 0;
    k = 0;
    while( cnt ~= r )
        k = k+1;
        if( s(ii(k),jj(k)) == -1 )
            cnt = cnt + 1;
        end
    end
    s(i,j) = 0;
    s( ii(k), jj(k) ) = 0;
    ri = ii(k);
    rj = jj(k);
    delta_pop = [-1, 2, -1];
    
% O site: CO + O --> CO2 + desorption
elseif( s(i,j) == -1 )
    suma = (s(i,ju)>0)  + (s(i,jd)>0)  + (s(il,j)>0)  + (s(ir,j)>0);
    
    r = ceil(suma*rand);
    cnt = 0;
    k = 0;
    while( cnt ~= r )
        k = k + 1;
        if( s(ii(k),jj(k)) == 1 )
            cnt = cnt+1;
        end
    end
    s(i,j) = 0;
    s( ii(k), jj(k) ) = 0; 
    ri = ii(k); 
    rj = jj(k);
    delta_pop = [-1, 2, -1];
else
    ri = ii(k); 
    rj = jj(k);
    delta_pop = [0, 0, 0];
end
    
delta_pop = delta_pop';
