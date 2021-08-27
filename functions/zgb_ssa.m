% File              : zgb_ssa.m
% Author            : George Arampatzis <garampat@ethz.ch>
% Date              : 27.08.2021
% Last Modified Date: 27.08.2021
% Last Modified By  : George Arampatzis <garampat@ethz.ch>
function [ts, s]  = zgb_ssa(data)
%% Descritpion:
%   Simulate the ZGB CO-oxidation model using the SSA algorithm.
%
%% Input:
%     data.s   :   initial lattice: square matrix with values {-1,0,+1}
%     data.rcd :   data struct containing information for the rates (see 'rates function')
%     data.t   :   sample the system at these times
% 
%% OUTPUT:
%     ts :   time series of coverage at jump times
%     s  :   lattice at final time

Nt = length(data.t);
T = data.t(end);

s = data.s;
if size(s,1) ~= size(s,2)
   error('The lattice must be square.')
end
N = size(s,1);

ts = zeros(3, Nt); 
ts(:,1) = observable(s);
ts_old = ts(:,1);

% compute initial rates on every site of the lattice
c = zeros(N,N);
for i=1:N
    for j=1:N
        c(i,j) = total_rate(s, i, j, data.rcd);
    end
end

ccum = cumsum(c(:));
c0 = ccum(end);

t = 0;
cnt_obs = 1;
delta_pop = zeros(3,1);

while t < T 
    
    % check for degenerate states
    if c0 <= 0 
        fprintf('\nSystem has reached an absorbing State. All rates are zero.\n');
        n = Nt - cnt_obs + 1;
        ts(:,cnt_obs:Nt) = repmat( ts_old, 1, n );
        return;
    end
    
    % find the site in lattice
    r = rand*c0 ;
    ind = find(ccum>r, 1, 'first');
    [ii, jj] = ind2sub( [N,N], ind);
    
    % compute waiting time
    dt = - log(rand)/c0;
    
    % store the observable before updating species
    ts_old =  ts_old + delta_pop;
    
    % move system to new state
    [ri, rj, s, delta_pop] = new_state(s, ii, jj, data.rcd);

    % update the observable array
    t_old = t ;
    t = t + dt;
    if( t_old <= data.t(cnt_obs) && t > data.t(cnt_obs)  )
        while( (cnt_obs <= Nt) && (data.t(cnt_obs) <= t) )
            ts(:,cnt_obs) = ts_old;
            cnt_obs = cnt_obs + 1;
        end
    end
    
    % compute new rates
    [il, ir, ju, jd] = neighbors_of(ii, jj, N);
    c(ii,jj) = total_rate( s, ii, jj, data.rcd );
    c(il,jj) = total_rate( s, il, jj, data.rcd );
    c(ir,jj) = total_rate( s, ir, jj, data.rcd );
    c(ii,ju) = total_rate( s, ii, ju, data.rcd );
    c(ii,jd) = total_rate( s, ii, jd, data.rcd );

    % if there was a reaction, update neighbors rates
    if(ri>0) 
        [il, ir, ju, jd] = neighbors_of( ri, rj, N );
        c(ri,rj) = total_rate( s, ri, rj, data.rcd );
        c(il,rj) = total_rate( s, il, rj, data.rcd );
        c(ir,rj) = total_rate( s, ir, rj, data.rcd );
        c(ri,ju) = total_rate( s, ri, ju, data.rcd );
        c(ri,jd) = total_rate( s, ri, jd, data.rcd ); 
    end
    
    % update cummulative sum of rates
    ccum = cumsum(c(:));
    c0 = ccum(end);
end

end

function ts = observable(s)
    ts = zeros( 3, 1 );
    ts(1) = sum( s(:)==-1 );
    ts(2) = sum( s(:)== 0 );
    ts(3) = sum( s(:)==+1 );
end
