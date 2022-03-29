clear;
setVerbLevs;
%setprngstates(39901520); % 39901520 oldCrude is particularly bad relative to fsovle.
%setprngstates(91565536); % 91565536 fsolve fails and is slow.
%setprngstates(0);
%setprngstates(95116336); %95116336 "advancing" is worse; probably a fluke.
%setprngstates(34444720); %34444720 Challenging 200x200
setprngstates();
numFigs = 0;
%
sizeX = 100;
sizeF = 100;
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
msg( __FILE__, __LINE__, sprintf( "rcond(matJE'*matJE) = %f.", rcond(matJE'*matJE) ) );
%
vecX0 = zeros(sizeX,1);
vecF0 = funchF(vecX0);
%
timeSS = time();
[ vecXF_125, vecFF_125, datOut_125 ] = findZero_125( vecX0, funchF );
time_125 = time()-timeSS;
%
timeSS = time();
[ vecXF_100, vecFF_100, datOut_100 ] = findZero_100( vecX0, funchF );
time_100 = time()-timeSS;
%
%
msg( __FILE__, __LINE__, sprintf( "norm(vecF0) = %g.", norm(vecF0) ) );
msg( __FILE__, __LINE__, "Solver results..." );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11s;  %11s;  %11s;  %11s.", "solver name", "||vecFF||", "iterCount", "fevalCount", "time(s)" ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "100", norm(vecFF_100), datOut_100.iterCount, datOut_100.fevalCount, time_100 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "125", norm(vecFF_125), datOut_125.iterCount, datOut_125.fevalCount, time_125 ) );
%
%
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_100.fevalCountVals, datOut_100.fNormVals+eps, '^-', 'markersize', 20, 'linewidth', 2, ...
  datOut_125.fevalCountVals, datOut_125.fNormVals+eps, 'v-', 'markersize', 20, 'linewidth', 2  );
grid on;
ylabel( "||f||" );
xlabel( "feval count" );
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_100.iterCountVals, datOut_100.fNormVals+eps, '^-', 'markersize', 20, 'linewidth', 2, ...
  datOut_125.iterCountVals, datOut_125.fNormVals+eps, 'v-', 'markersize', 20, 'linewidth', 2 );
grid on;
ylabel( "||f||" );
xlabel( "iteration count" );
