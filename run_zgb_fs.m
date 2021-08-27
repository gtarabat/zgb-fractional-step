% File              : run_zgb_fs.m
% Author            : George Arampatzis <garampat@ethz.ch>
% Date              : 27.08.2021
% Last Modified Date: 27.08.2021
% Last Modified By  : George Arampatzis <garampat@ethz.ch>
%
% Description: 
%   Run and plot the ZGB model using the Fractional-Step algorithm.

clear;clc

addpath('functions')

N = 100;
data.rcd.k1 = 0.394;
data.rcd.k2 = 1;
data.s = 0*round(rand(N,N));
data.nblocks = 4;
data.T = 10;
data.dt = 0.01;

video.save = true;
video.Nf = 10;

[t, ts, s, total_workload] = zgb_fs( data, video );


N2 = N^2;
figure();
plot( t, ts(1,:)/N2, '-'); hold on
plot( t, ts(2,:)/N2, '-')
plot( t, ts(3,:)/N2, '-')
grid on
axis tight

ax = gca;
ax.XLabel.String = 'time';
ax.YLabel.String = 'concentration';

legend('$O$','$\emptyset$','$C)$', 'Interpreter','latex')

rmpath('functions')
