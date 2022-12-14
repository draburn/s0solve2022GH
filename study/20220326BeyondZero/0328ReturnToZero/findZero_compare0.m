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
[ vecXF_200, vecFF_200, datOut_200 ] = findZero_200( vecX0, funchF );
time_200 = time()-timeSS;
%return
%
timeSS = time();
[ vecXF_fsolve, vecFF_fsolve, datOut_fsolve, ] = findZero_fsolve( vecX0, funchF );
time_fsolve = time()-timeSS;
[ vecXF_fg, vecFF_fg, datOut_fg ] = findZero_fsolveGnostic( vecX0, funchF );
%
timeSS = time();
[ vecXF_000, vecFF_000, datOut_000 ] = findZero_000( vecX0, funchF );
time_000 = time()-timeSS;
%
timeSS = time();
[ vecXF_050, vecFF_050, datOut_050 ] = findZero_050( vecX0, funchF );
time_050 = time()-timeSS;
%
timeSS = time();
[ vecXF_090, vecFF_090, datOut_090 ] = findZero_090( vecX0, funchF );
time_090 = time()-timeSS;
%
timeSS = time();
[ vecXF_100, vecFF_100, datOut_100 ] = findZero_100( vecX0, funchF );
time_100 = time()-timeSS;
%
%
msg( __FILE__, __LINE__, sprintf( "norm(vecF0) = %g.", norm(vecF0) ) );
msg( __FILE__, __LINE__, "Solver results..." );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11s;  %11s;  %11s;  %11s.", "solver name", "||vecFF||", "iterCount", "fevalCount", "time(s)" ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "fsolve", norm(vecFF_fsolve), datOut_fsolve.iterCount, datOut_fsolve.fevalCount, time_fsolve ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "000", norm(vecFF_000), datOut_000.iterCount, datOut_000.fevalCount, time_000 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "050", norm(vecFF_050), datOut_050.iterCount, datOut_050.fevalCount, time_050 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "090", norm(vecFF_090), datOut_090.iterCount, datOut_090.fevalCount, time_090 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "100", norm(vecFF_100), datOut_100.iterCount, datOut_100.fevalCount, time_100 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "200", norm(vecFF_200), datOut_200.iterCount, datOut_200.fevalCount, time_200 ) );
%
%
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_fg.fevalCountVals, datOut_fg.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_000.fevalCountVals, datOut_000.fNormVals+eps, '+-', 'markersize', 20, 'linewidth', 2, ...
  datOut_050.fevalCountVals, datOut_050.fNormVals+eps, 'x-', 'markersize', 20, 'linewidth', 2, ...
  datOut_090.fevalCountVals, datOut_090.fNormVals+eps, 's-', 'markersize', 20, 'linewidth', 2, ...
  datOut_100.fevalCountVals, datOut_100.fNormVals+eps, '^-', 'markersize', 20, 'linewidth', 2, ...
  datOut_200.fevalCountVals, datOut_200.fNormVals+eps, 'v-', 'markersize', 20, 'linewidth', 2  );
grid on;
ylabel( "||f||" );
xlabel( "feval count" );
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_fg.iterCountVals, datOut_fg.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_000.iterCountVals, datOut_000.fNormVals+eps, '+-', 'markersize', 20, 'linewidth', 2, ...
  datOut_050.iterCountVals, datOut_050.fNormVals+eps, 'x-', 'markersize', 20, 'linewidth', 2, ...
  datOut_090.iterCountVals, datOut_090.fNormVals+eps, 's-', 'markersize', 20, 'linewidth', 2, ...
  datOut_100.iterCountVals, datOut_100.fNormVals+eps, '^-', 'markersize', 20, 'linewidth', 2, ...
  datOut_200.iterCountVals, datOut_200.fNormVals+eps, 'v-', 'markersize', 20, 'linewidth', 2 );
grid on;
ylabel( "||f||" );
xlabel( "iteration count" );
