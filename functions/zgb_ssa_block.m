% File              : zgb_ssa_block.m
% Author            : George Arampatzis <garampat@ethz.ch>
% Date              : 27.08.2021
% Last Modified Date: 27.08.2021
% Last Modified By  : George Arampatzis <garampat@ethz.ch>
function [s, ts, events]  = zgb_ssa_block(s, T, rcd, is, js, n)
%% Descritpion:
%   Simulate the ZGB CO-oxidation model using the SSA algorithm on a
%   sublattice of s. The sublattice starts at the site (is,js) and has
%   length n on the right and length n on the down.
%
%% Input:
%     s   :   initial lattice: square matrix with values {-1,0,+1}
%     T   :   simulate up to time T
%     rcd :   data struct containing information for the rates (see 'rates function')
%     is  :   starting x-position of the sub-lattice
%     js  :   starting y-position of the sub-lattice
%     n   :   side length of the sub-lattice   
%
%% OUTPUT:
%     s  :   lattice at final time
%     ts :   time series of coverage at jump times
%     events : numbers of reactions took place until time T


c = zeros(n, n);

for i = 1:n
    for j = 1:n
        c(i,j) = total_rate(s, is+i, js+j, rcd);
    end
end

tmp = c';
ccum = cumsum(tmp(:));
c0 = ccum(end);

ts = zeros(3,1); 
ts(1) = sum(s(is+1:is+n, js+1:js+n) == -1, 'all');
ts(2) = sum(s(is+1:is+n, js+1:js+n) ==  0, 'all');
ts(3) = sum(s(is+1:is+n, js+1:js+n) == +1, 'all');

t = -log(rand) / c0;

tt = t;
events = 1;

while tt <= T 
    
    % find the site in the local lattice
    r = rand * c0;
    ind = find(ccum>r, 1, 'first');
    [jj, ii] = ind2sub( [n,n], ind);
        
    % indices of the site in the global lattice
    sii = is + ii; 
    sjj = js + jj;
    
    % move system to new state
    [sri, srj, s, popCh] = new_state(s, sii, sjj, rcd);
    
    ts = ts + popCh;
    
    % compute new total_rate only on the sites of the local block. If the
    % neighbours belong in neighbouring blocks, do not update their total_rate
    c(ii,jj) = total_rate(s, sii, sjj, rcd);
    if jj-1 > 0
        c(ii,jj-1) = total_rate(s, sii, sjj-1, rcd);
    end
    if jj+1 <= n
        c(ii,jj+1) = total_rate(s, sii, sjj+1, rcd);
    end
    if ii-1 > 0
        c(ii-1,jj) = total_rate(s, sii-1, sjj, rcd);
    end
    if ii+1 <= n
        c(ii+1,jj) = total_rate(s, sii+1, sjj, rcd);
    end
    
    % if there was a reaction, update neighbors total_rate
    ri = sri - is; 
    rj = srj - js;
    
    if( sri > 0 && srj > 0 )
    
        if( ri > 0 && ri <= n && rj > 0 && rj <= n )
            c(ri,rj) = total_rate(s, sri, srj, rcd);
        end
        if( ri > 0 && ri <= n && rj-1 > 0 && rj-1 <= n )
            c(ri,rj-1) = total_rate(s, sri, srj-1, rcd);
        end
        if( ri > 0 && ri <= n && rj+1 > 0 && rj+1 <= n )
            c(ri,rj+1) = total_rate(s,sri,srj+1, rcd);
        end
        if( ri-1 > 0 && ri-1 <= n && rj > 0 && rj <= n )
            c(ri-1,rj) = total_rate(s, sri-1, srj, rcd);
        end
        if( ri+1 > 0 && ri+1 <= n && rj > 0 && rj <= n )
            c(ri+1,rj) = total_rate(s, sri+1, srj, rcd);
        end

    end
    
    tmp = c';
    ccum = cumsum(tmp(:));
    c0 = ccum(end);
    
    t = -log(rand) / c0;
    
    tt = tt + t;
    events = events + 1;
end

    
