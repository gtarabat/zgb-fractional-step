% File              : run_zgb_fs.m
% Author            : George Arampatzis <garampat@ethz.ch>
% Date              : 27.08.2021
% Last Modified Date: 27.08.2021
% Last Modified By  : George Arampatzis <garampat@ethz.ch>
%
% Description:
%   Compare the Fractional-Step algorithm with SSA on the ZGB model

clear;clc

addpath('functions')

N = 10;
T = 8;
s0 = zeros(N, N);
rcd.k1 = 0.394;
rcd.k2 = 1;

N2 = N^2;

%% Run Fractional-Step

data.rcd = rcd;
data.s = s0;
data.nblocks = 2;
data.T = T;
data.dt = 0.01;

video.save = false;

[t, ts, s, total_workload] = zgb_fs( data, video );

figure();
plot( t, ts(1,:)/N2, '-'); hold on
plot( t, ts(2,:)/N2, '-')
plot( t, ts(3,:)/N2, '-')
grid on
axis tight

ax = gca;
ax.XLabel.String = 'time';
ax.YLabel.String = 'concentration';

%% Run SSA
data = [];

data.rcd = rcd;
data.s = s0;
data.t = 0: 0.01: T;

[ts, s] = zgb_ssa(data);

set(gca,'ColorOrderIndex',1)
plot( data.t, ts(1,:)/N2, '--' )
plot( data.t, ts(2,:)/N2, '--' )
plot( data.t, ts(3,:)/N2, '--' )
grid on
axis tight

legend('$O$ (FS)','$\emptyset$ (FS)','$C)$ (FS)',...
       '$O (SSA)$','$\emptyset$ (SSA)','$CO$ (SSA)', 'Interpreter','latex')

rmpath('functions')
