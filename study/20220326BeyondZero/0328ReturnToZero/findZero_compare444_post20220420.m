%%
%restartIter = 1 % Converges.
%restartIter = 5 % Cnvg
%restartIter = 8 % Cnvg
%%restartIter = 9 % Stall
%restartIter = 10 % Stalls.
%restartIter = 11 % Stall
restartIter = 15 % Stall
%restartIter = 16 %Stall
%restartIter = 17 %Stall
%restartIter = 18 % Cnvgish
%restartIter = 19 % Cnvgish
%restartIter = 20 % Converge slowly.
%restartIter = 21 % Converge slowly.
%restartIter = 50 % Converges-ish slowly
%%
%vecX_restart = datOut_800.vecXVals(:,restartIter);
alpha = 0.5
vecX_restart = (1.0+alpha)*datOut_800.vecXVals(:,restartIter) - alpha*vecX0;
%
timeSS = time();
prm = [];
prm.iterMax = 200;
[ vecXF_post, vecFF_post, datOut_post ] = findZero_050( vecX_restart, funchF, prm );
time_post = time()-timeSS
%
%
%
numFigs = 220;
numFigs++; figure( numFigs );
semilogy( ...
  datOut_800.fevalCountVals, datOut_800.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_post.fevalCountVals, datOut_post.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_550.fevalCountVals, datOut_550.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_444.fevalCountVals, datOut_444.fNormVals+eps, '^-', 'markersize', 20, 'linewidth', 2    );
grid on;
ylabel( "||f||" );
xlabel( "feval count" );
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_800.iterCountVals, datOut_800.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_post.iterCountVals, datOut_post.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_550.iterCountVals, datOut_550.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_444.iterCountVals, datOut_444.fNormVals+eps, '^-', 'markersize', 20, 'linewidth', 2 );
grid on;
ylabel( "||f||" );
xlabel( "iteration count" );
