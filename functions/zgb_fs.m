% File              : zgb_fs.m
% Author            : George Arampatzis <garampat@ethz.ch>
% Date              : 27.08.2021
% Last Modified Date: 27.08.2021
% Last Modified By  : George Arampatzis <garampat@ethz.ch>
function [t, ts, s, total_workload] = zgb_fs( data, video)
%% Descritpion:
%   Simulate the ZGB CO-oxidation model using the Fractional-Step algorithm.
%   The implementation is not parallel and it runs all blocks sequantially.
%
%% Input:
%     data.s   :   initial lattice: square matrix with values {-1,0,+1}
%     data.rcd :   data struct containing information for the rates (see 'rates function')
%     data.T   :   simulate the system up to this time
%     data.dt  :   time step for the Fractional-Step algorithm
%     data.nblocks : number of blocks on each direction, must be even
%
%     video.save :  boolean variable to save a video of the simulation
%     video.Nf   :  save the lattice in the video every Nf steps
% 
%% OUTPUT:
%     t  :   vector with time instances
%     ts :   time series of coverage at jump times
%     s  :   lattice at final time
%     total_workload: the total number of events for every block

s = data.s;

if size(s,1) ~= size(s,2)
   error('The lattice must be square.') 
end

N = size(s,1);
nblocks = data.nblocks;

if mod(nblocks,2) ~= 0
    error('The number of blocks must be even.');
end

if mod(N, nblocks) ~= 0
    error('The length of the laticce side must be devided by the number of blocks.');
end

if video.save == true
    vw = VideoWriter('zgb_fs.avi');
    open(vw);
    fig = figure(1);
    fig.Position = [636, 24, 1562, 1313];
    colormap( [1 0 0; 1 1 1; 0 0 0] );
    ax1 = subplot(5,1,1:4);
    ax1.Visible = 'off';
    set( findall( ax1, 'type', 'text'), 'visible', 'on')
    ax2 = subplot(5,1,5);
end

Ns = data.T/data.dt + 1;
sblock = N/nblocks;

workload = zeros(nblocks,nblocks);
total_workload = workload;

JIND = zeros(2,nblocks/2,2);
JIND(:,:,1) = [1:2:nblocks ; 2:2:nblocks];
JIND(:,:,2) = [2:2:nblocks ; 1:2:nblocks];

ts = zeros(3,Ns); 
ts(1,1) = sum(s(:)==-1);
ts(2,1) = sum(s(:)== 0);
ts(3,1) = sum(s(:)==+1);

for step = 2:Ns
   
    for k = 1:2
    for i = 1:nblocks
    for j = JIND(2-mod(i,2),:,k)
        
        is = (i-1)*sblock;
        js = (j-1)*sblock;
        
        [s, ts_local, workload(i,j)]  = zgb_ssa_block(s, data.dt, data.rcd, is, js, sblock);
        
        total_workload(i,j) = total_workload(i,j) + workload(i,j);
        ts(:,step) = ts(:,step) + ts_local;
    end
    end
    end
    
    % add the lattice to the video
    if( video.save == true && mod(step,video.Nf) == 0 )
        cla(ax1)
        imagesc(ax1, s); 
        hold(ax1,'on')
        ax1.Visible = 'off';
        set( findall( ax1, 'type', 'text'), 'visible', 'on')
        ax1.Title.String = [ 't = ' num2str(step*data.dt) ];
        for i=1:nblocks
            for j=1:nblocks
                plot(ax1, [(j-1)*sblock j*sblock]+0.5, [  i  *sblock i*sblock]+0.5, 'b.-');
                plot(ax1, [  j  *sblock j*sblock]+0.5, [(i-1)*sblock i*sblock]+0.5, 'b.-');
            end
        end
        axis(ax1, 'equal')
        
        wl = reshape(workload',1,nblocks^2);
        
        cla(ax2)
        bar(ax2, 1:nblocks^2, wl); 
        xlim([0, nblocks^2+1])
        ax2.Title.String = 'Workload';
        ax2.XLabel.String = '# block';
        ax2.YLabel.String = 'number of reactions per block';
        
        frame = getframe(fig);
        writeVideo(vw, frame);
    end
    
end

t = (0:Ns-1)*data.dt;

if video.save == true
    close(vw);
end
