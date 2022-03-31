clear;
setVerbLevs;
%setprngstates(0);
%setprngstates(90186240); %90186240 for 150x150. Some stall.
%setprngstates(76252688); %76252688 for 150x150 causes div 0 in hack code in lisolf_directed.
%setprngstates(53244560); % 550 very slow; 880 okay.
setprngstates(0); % 800 stalls but 550 does not.
%setprngstates(9216528); % 800 stalls but 550 does not.
%setprngstates(17374144); % 800 is weird but converges.
%setprngstates(53858336); % 800 stalls but 550 does not.
numFigs = 0;
%
sizeX = 150; sizeF = 150;
%
vecXE = randn(sizeX,1);
matJE = diag(1.0+abs(randn(min([sizeF,sizeX]),1))) + 0.3*randn(sizeF,sizeX);
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
%vecX0 = zeros(sizeX,1);
vecX0 = vecXE + 0.5*randn(sizeX,1);
vecF0 = funchF(vecX0);
Df = jacobs( vecX0, funchF );
msg( __FILE__, __LINE__, sprintf( "rcond(Df'*Df) = %0.3e.", rcond(Df'*Df) ) );
%
prm = [];
timeSS = time();
[ vecXF_444, vecFF_444, datOut_444 ] = findZero_444( vecX0, funchF, prm );
time_444 = time()-timeSS;
%
%
prm = [];
prm.iterMax = 5;
prm.iterMax = 200;
timeSS = time();
[ vecXF_550, vecFF_550, datOut_550 ] = findZero_550( vecX0, funchF, prm );
time_550 = time()-timeSS;
%
%
prm = [];
prm.iterMax = 200;
timeSS = time();
[ vecXF_800, vecFF_800, datOut_800 ] = findZero_800( vecX0, funchF, prm );
time_800 = time()-timeSS;
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
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "800", norm(vecFF_800), datOut_800.iterCount, datOut_800.fevalCount, time_800 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "550", norm(vecFF_550), datOut_550.iterCount, datOut_550.fevalCount, time_550 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "444", norm(vecFF_444), datOut_444.iterCount, datOut_444.fevalCount, time_444 ) );
%
%
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_800.fevalCountVals, datOut_800.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_550.fevalCountVals, datOut_550.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_444.fevalCountVals, datOut_444.fNormVals+eps, '^-', 'markersize', 20, 'linewidth', 2    );
grid on;
ylabel( "||f||" );
xlabel( "feval count" );
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_800.iterCountVals, datOut_800.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_550.iterCountVals, datOut_550.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_444.iterCountVals, datOut_444.fNormVals+eps, '^-', 'markersize', 20, 'linewidth', 2 );
grid on;
ylabel( "||f||" );
xlabel( "iteration count" );
