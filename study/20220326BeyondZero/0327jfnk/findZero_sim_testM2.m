clear;
setVerbLevs;
setprngstates(0);
numFigs = 0;
%
sizeX = 20;
sizeF = 20;
funchF = @(x)( x - 1.0 );
%
vecX0 = zeros(sizeX,1);
%
timeSS = time();
[ vecXF_fsolve, vecFF_fsolve, datOut_fsolve, ] = findZero_fsolve( vecX0, funchF );
time_fsolve = time()-timeSS;
[ vecXF_fg, vecFF_fg, datOut_fg ] = findZero_fsolveGnostic( vecX0, funchF );
%
timeSS = time();
[ vecXF_sim, vecFF_sim, datOut_sim ] = findZero_shouldIterMin( vecX0, funchF );
time_sim = time()-timeSS;
%
vecF0 = funchF(vecX0);
msg( __FILE__, __LINE__, sprintf( "norm(vecF0) = %g.", norm(vecF0) ) );
msg( __FILE__, __LINE__, "Solver results..." );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11s;  %11s;  %11s;  %11s.", "solver name", "||vecFF||", "iterCount", "fevalCount", "time(s)" ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "fsolve", norm(vecFF_fsolve), datOut_fsolve.iterCount, datOut_fsolve.fevalCount, time_fsolve ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "shouldIterMin", norm(vecFF_sim), datOut_sim.iterCount, datOut_sim.fevalCount, time_sim ) );
%
%
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_fg.fevalCountVals, datOut_fg.fNormVals+eps^2, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_sim.fevalCountVals, datOut_sim.fNormVals+eps^2, '+-', 'markersize', 20, 'linewidth', 2 );
grid on;
ylabel( "||f||" );
xlabel( "feval count" );
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_fg.iterCountVals, datOut_fg.fNormVals+eps^2, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_sim.iterCountVals, datOut_sim.fNormVals+eps^2, '+-', 'markersize', 20, 'linewidth', 2 );
grid on;
ylabel( "||f||" );
xlabel( "iteration count" );
