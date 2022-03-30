clear;
setVerbLevs;
%setprngstates(0);
%setprngstates(90186240); %90186240 for 150x150. Some stall.
setprngstates(76252688); %76252688 for 150x150 causes div 0 in hack code in lisolf_directed.
numFigs = 0;
%
sizeX = 150; sizeF = 150;
%sizeX = 15; sizeF = 15;
%sizeX = 5; sizeF = 5;
%
vecXE = randn(sizeX,1);
matJE = diag(1.0+abs(randn(min([sizeF,sizeX]),1))) + 0.01*randn(sizeF,sizeX);
matA0 = 0.0001*randn(sizeF,sizeX);
matA1 = randn(sizeX,sizeX);
matA2 = randn(sizeX,sizeX);
matB0 = 0.0001*randn(sizeF,sizeX);
matB1 = randn(sizeX,sizeX);
matB2 = randn(sizeX,sizeX);
matB3 = randn(sizeX,sizeX);
y = @(x)( x - vecXE );
funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
msg( __FILE__, __LINE__, sprintf( "rcond(matJE'*matJE) = %0.3e.", rcond(matJE'*matJE) ) );
%
vecX0 = zeros(sizeX,1);
vecF0 = funchF(vecX0);
Df = jacobs( vecX0, funchF );
msg( __FILE__, __LINE__, sprintf( "rcond(Df'*Df) = %0.3e.", rcond(Df'*Df) ) );
%
prm.iterMax = 200;
%
timeSS = time();
[ vecXF_825, vecFF_825, datOut_825 ] = findZero_825( vecX0, funchF, prm );
time_825 = time()-timeSS;
%
timeSS = time();
[ vecXF_800, vecFF_800, datOut_800 ] = findZero_800( vecX0, funchF, prm );
time_800 = time()-timeSS;
%
timeSS = time();
[ vecXF_700, vecFF_700, datOut_700 ] = findZero_700( vecX0, funchF, prm );
time_700 = time()-timeSS;
%
timeSS = time();
[ vecXF_550, vecFF_550, datOut_550 ] = findZero_550( vecX0, funchF, prm );
time_550 = time()-timeSS;
%
timeSS = time();
[ vecXF_200, vecFF_200, datOut_200 ] = findZero_200( vecX0, funchF, prm );
time_200 = time()-timeSS;
%
timeSS = time();
[ vecXF_fsolve, vecFF_fsolve, datOut_fsolve, ] = findZero_fsolve( vecX0, funchF );
time_fsolve = time()-timeSS;
%
%
msg( __FILE__, __LINE__, sprintf( "norm(vecF0) = %g.", norm(vecF0) ) );
msg( __FILE__, __LINE__, "Solver results..." );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11s;  %11s;  %11s;  %11s.", "solver name", "||vecFF||", "iterCount", "fevalCount", "time(s)" ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "fsolve", norm(vecFF_fsolve), datOut_fsolve.iterCount, datOut_fsolve.fevalCount, time_fsolve ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "200", norm(vecFF_200), datOut_200.iterCount, datOut_200.fevalCount, time_200 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "550", norm(vecFF_550), datOut_550.iterCount, datOut_550.fevalCount, time_550 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "700", norm(vecFF_700), datOut_700.iterCount, datOut_700.fevalCount, time_700 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "800", norm(vecFF_800), datOut_800.iterCount, datOut_800.fevalCount, time_800 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "825", norm(vecFF_825), datOut_825.iterCount, datOut_825.fevalCount, time_825 ) );
%
%
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_200.fevalCountVals, datOut_200.fNormVals+eps, 'v-', 'markersize', 20, 'linewidth', 2, ...
  datOut_550.fevalCountVals, datOut_550.fNormVals+eps, '+-', 'markersize', 20, 'linewidth', 2, ...
  datOut_700.fevalCountVals, datOut_700.fNormVals+eps, 's-', 'markersize', 20, 'linewidth', 2, ...
  datOut_800.fevalCountVals, datOut_800.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_825.fevalCountVals, datOut_825.fNormVals+eps, '^-', 'markersize', 20, 'linewidth', 2    );
grid on;
ylabel( "||f||" );
xlabel( "feval count" );
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_200.iterCountVals, datOut_200.fNormVals+eps, 'v-', 'markersize', 20, 'linewidth', 2, ...
  datOut_550.iterCountVals, datOut_550.fNormVals+eps, '+-', 'markersize', 20, 'linewidth', 2, ...
  datOut_700.iterCountVals, datOut_700.fNormVals+eps, 's-', 'markersize', 20, 'linewidth', 2, ...
  datOut_800.iterCountVals, datOut_800.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_825.iterCountVals, datOut_825.fNormVals+eps, '^-', 'markersize', 20, 'linewidth', 2 );
grid on;
ylabel( "||f||" );
xlabel( "iteration count" );
