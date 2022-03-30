clear;
setVerbLevs;
setprngstates(0);
%setprngstates(39901520); % 39901520 oldCrude is particularly bad relative to fsovle.
%setprngstates(91565536); % 91565536 fsolve fails and is slow.
%setprngstates(95116336); %95116336 "advancing" is worse; probably a fluke.
%setprngstates(34444720); %34444720 Challenging 200x200
%setprngstates(8226384); %8226384 triggers a path.
%setprngstates();
numFigs = 0;
%
sizeX = 100;
sizeF = 100;
%
vecXE = randn(sizeX,1);
matJE = diag(1.0+abs(randn(min([sizeF,sizeX]),1))) + 0.001*randn(sizeF,sizeX);
%%matA0 = 0.00001*randn(sizeF,sizeX);
matA0 = 0.0001*randn(sizeF,sizeX);
matA1 = randn(sizeX,sizeX);
matA2 = randn(sizeX,sizeX);
%%matB0 = 0.00001*randn(sizeF,sizeX);
matB0 = 0.0001*randn(sizeF,sizeX);
matB1 = randn(sizeX,sizeX);
matB2 = randn(sizeX,sizeX);
matB3 = randn(sizeX,sizeX);
y = @(x)( x - vecXE );
funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
msg( __FILE__, __LINE__, sprintf( "rcond(matJE'*matJE) = %f.", rcond(matJE'*matJE) ) );
%
%vecX0 = zeros(sizeX,1);
vecX0 = vecXE + 0.8*randn(sizeX,1);
vecF0 = funchF(vecX0);
prm.iterMax = 200;
%
timeSS = time();
[ vecXF_105, vecFF_105, datOut_105 ] = findZero_105( vecX0, funchF, prm );
time_105 = time()-timeSS;
%
timeSS = time();
[ vecXF_700, vecFF_700, datOut_700 ] = findZero_700( vecX0, funchF, prm );
time_700 = time()-timeSS;
%
timeSS = time();
[ vecXF_600, vecFF_600, datOut_600 ] = findZero_600( vecX0, funchF, prm );
time_600 = time()-timeSS;
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
[ vecXF_100, vecFF_100, datOut_100 ] = findZero_100( vecX0, funchF, prm );
time_100 = time()-timeSS;
%
timeSS = time();
[ vecXF_fsolve, vecFF_fsolve, datOut_fsolve, ] = findZero_fsolve( vecX0, funchF );
time_fsolve = time()-timeSS;
[ vecXF_fg, vecFF_fg, datOut_fg ] = findZero_fsolveGnostic( vecX0, funchF );
%
%
msg( __FILE__, __LINE__, sprintf( "norm(vecF0) = %g.", norm(vecF0) ) );
msg( __FILE__, __LINE__, "Solver results..." );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11s;  %11s;  %11s;  %11s.", "solver name", "||vecFF||", "iterCount", "fevalCount", "time(s)" ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "fsolve", norm(vecFF_fsolve), datOut_fsolve.iterCount, datOut_fsolve.fevalCount, time_fsolve ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "100", norm(vecFF_100), datOut_100.iterCount, datOut_100.fevalCount, time_100 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "200", norm(vecFF_200), datOut_200.iterCount, datOut_200.fevalCount, time_200 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "550", norm(vecFF_550), datOut_550.iterCount, datOut_550.fevalCount, time_550 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "600", norm(vecFF_600), datOut_600.iterCount, datOut_600.fevalCount, time_600 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "700", norm(vecFF_700), datOut_700.iterCount, datOut_700.fevalCount, time_700 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "105", norm(vecFF_105), datOut_105.iterCount, datOut_105.fevalCount, time_105 ) );
%
%
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_fg.fevalCountVals, datOut_fg.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_100.fevalCountVals, datOut_100.fNormVals+eps, '^-', 'markersize', 20, 'linewidth', 2, ...
  datOut_200.fevalCountVals, datOut_200.fNormVals+eps, 'v-', 'markersize', 20, 'linewidth', 2, ...
  datOut_550.fevalCountVals, datOut_550.fNormVals+eps, '+-', 'markersize', 20, 'linewidth', 2, ...
  datOut_600.fevalCountVals, datOut_600.fNormVals+eps, 'x-', 'markersize', 20, 'linewidth', 2, ...
  datOut_700.fevalCountVals, datOut_700.fNormVals+eps, 's-', 'markersize', 20, 'linewidth', 2, ...
  datOut_105.fevalCountVals, datOut_105.fNormVals+eps, 's-', 'markersize', 20, 'linewidth', 2   );
grid on;
ylabel( "||f||" );
xlabel( "feval count" );
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_fg.iterCountVals, datOut_fg.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_100.iterCountVals, datOut_100.fNormVals+eps, '^-', 'markersize', 20, 'linewidth', 2, ...
  datOut_200.iterCountVals, datOut_200.fNormVals+eps, 'v-', 'markersize', 20, 'linewidth', 2, ...
  datOut_550.iterCountVals, datOut_550.fNormVals+eps, '+-', 'markersize', 20, 'linewidth', 2, ...
  datOut_600.iterCountVals, datOut_600.fNormVals+eps, 'x-', 'markersize', 20, 'linewidth', 2, ...
  datOut_700.iterCountVals, datOut_700.fNormVals+eps, 's-', 'markersize', 20, 'linewidth', 2, ...
  datOut_105.iterCountVals, datOut_105.fNormVals+eps, 's-', 'markersize', 20, 'linewidth', 2 );
grid on;
ylabel( "||f||" );
xlabel( "iteration count" );
