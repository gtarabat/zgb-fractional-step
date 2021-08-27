% File              : run_zgb_ssa.m
% Author            : George Arampatzis <garampat@ethz.ch>
% Date              : 27.08.2021
% Last Modified Date: 27.08.2021
% Last Modified By  : George Arampatzis <garampat@ethz.ch>
%
% Description: 
%   Run and plot the ZGB model using the SSA algorithm

clear;clc

addpath('functions')

N = 100;
data.rcd.k1 = 0.394;
data.rcd.k2 = 1;
data.s = 0*round(rand(N,N));
data.t = 0: 0.01: 100;

[ts, s ] = zgb_ssa(data);

N2 = N^2;
figure()
hold on
plot( data.t, ts(1,:)/N2, '-' )
plot( data.t, ts(2,:)/N2, '-' )
plot( data.t, ts(3,:)/N2, '-' )
grid on
axis tight

legend('$O$','$\emptyset$','$CO$', 'Interpreter','latex')

rmpath('functions')
