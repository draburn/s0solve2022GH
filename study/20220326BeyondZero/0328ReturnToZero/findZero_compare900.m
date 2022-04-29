clear;
setVerbLevs;
%setprngstates(0);
%setprngstates(90186240); %90186240 for 150x150. Some stall.
%setprngstates(76252688); %76252688 for 150x150 causes div 0 in hack code in lisolf_directed.
%setprngstates(53244560); % 550 very slow; 880 okay.
%setprngstates(0); % 800 stalls but 550 does not.
%setprngstates(9216528); % 800 stalls but 550 does not.
%setprngstates(17374144); % 800 is weird but converges.
%setprngstates(53858336); % 800 stalls but 550 does not.
%setprngstates(48689984); % 904s are bad, others are good; note: even fsolve stalls.
% Two past 48689984?: almost everything stalls.
%setprngstates(36693936); % Most stall, but not fzero.
setprngstates(22931408); % 940x200 is bad for 150x150.
%setprngstates();
numFigs = 0;
%
sizeX = 150; sizeF = 150;
%%%sizeX = 100; sizeF = 100;
%%%sizeX = 50; sizeF = 50;
%%%sizeX = 20; sizeF = 20;
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
timeSS = time();
[ vecXF_fsolve, vecFF_fsolve, datOut_fsolve, ] = findZero_fsolve( vecX0, funchF );
time_fsolve = time()-timeSS;
msg( __FILE__, __LINE__, "fsolve results..." );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "fsolve", norm(vecFF_fsolve), datOut_fsolve.iterCount, datOut_fsolve.fevalCount, time_fsolve ) );
%
prm = [];
prm.iterMax = 50;
prm.slinsolfver = 200;
timeSS = time();
[ vecXF_940x200, vecFF_940x200, datOut_940x200 ] = findZero_940( vecX0, funchF, prm );
time_940x200 = time()-timeSS;
%
%
prm = [];
prm.iterMax = 10;
prm.slinsolfver = 100;
timeSS = time();
[ vecXF_940, vecFF_940, datOut_940 ] = findZero_940( vecX0, funchF, prm );
time_940 = time()-timeSS;
%
%
%
prm = [];
prm.step_prm.usePostLinsolfPhiPatch = false;
prm.iterMax = 5;
prm.iterMax = 200;
timeSS = time();
[ vecXF_904, vecFF_904, datOut_904 ] = findZero_904( vecX0, funchF, prm );
time_904 = time()-timeSS;
%
%
prm = [];
prm.step_prm.usePostLinsolfPhiPatch = true;
prm.iterMax = 5;
prm.iterMax = 200;
timeSS = time();
[ vecXF_904x, vecFF_904x, datOut_904x ] = findZero_904( vecX0, funchF, prm );
time_904x = time()-timeSS;
%
prm = [];
timeSS = time();
[ vecXF_444, vecFF_444, datOut_444 ] = findZero_444( vecX0, funchF, prm );
time_444 = time()-timeSS;
%
prm = [];
prm.useCoasting = false;
timeSS = time();
[ vecXF_444x, vecFF_444x, datOut_444x ] = findZero_444( vecX0, funchF, prm );
time_444x = time()-timeSS;
%
%
prm = [];
prm.iterMax = 200;
prm.linsolf_tol0 = 0.01;
timeSS = time();
[ vecXF_900, vecFF_900, datOut_900 ] = findZero_900( vecX0, funchF, prm );
time_900 = time()-timeSS;
%
%
prm = [];
prm.iterMax = 200;
timeSS = time();
[ vecXF_800, vecFF_800, datOut_800 ] = findZero_800( vecX0, funchF, prm );
time_800 = time()-timeSS;
%
%










% Push to triple-digit line number.
msg( __FILE__, __LINE__, sprintf( "norm(vecF0) = %g.", norm(vecF0) ) );
msg( __FILE__, __LINE__, "Solver results..." );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11s;  %11s;  %11s;  %11s.", "solver name", "||vecFF||", "iterCount", "fevalCount", "time(s)" ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "fsolve", norm(vecFF_fsolve), datOut_fsolve.iterCount, datOut_fsolve.fevalCount, time_fsolve ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "800", norm(vecFF_800), datOut_800.iterCount, datOut_800.fevalCount, time_800 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "900", norm(vecFF_900), datOut_900.iterCount, datOut_900.fevalCount, time_900 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "940", norm(vecFF_940), datOut_940.iterCount, datOut_940.fevalCount, time_940 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "940x200", norm(vecFF_940x200), datOut_940x200.iterCount, datOut_940x200.fevalCount, time_940x200 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "904", norm(vecFF_904), datOut_904.iterCount, datOut_904.fevalCount, time_904 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "904x", norm(vecFF_904x), datOut_904x.iterCount, datOut_904x.fevalCount, time_904x ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "444", norm(vecFF_444), datOut_444.iterCount, datOut_444.fevalCount, time_444 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "444x", norm(vecFF_444x), datOut_444x.iterCount, datOut_444x.fevalCount, time_444x ) );
%
%
%
epsViz = 1.0e-18;
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_800.fevalCountVals, datOut_800.fNormVals+epsViz, 'x-', 'markersize', 20, 'linewidth', 2, ...
  datOut_900.fevalCountVals, datOut_900.fNormVals+epsViz, '+-', 'markersize', 20, 'linewidth', 2, ...
  datOut_940.fevalCountVals, datOut_940.fNormVals+epsViz, '*-', 'markersize', 20, 'linewidth', 2, ...
  datOut_940x200.fevalCountVals, datOut_940x200.fNormVals+epsViz, 'p-', 'markersize', 20, 'linewidth', 2, ...
  datOut_904.fevalCountVals, datOut_904.fNormVals+epsViz, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_904x.fevalCountVals, datOut_904x.fNormVals+epsViz, 's-', 'markersize', 20, 'linewidth', 2, ...
  datOut_444.fevalCountVals, datOut_444.fNormVals+epsViz, '^-', 'markersize', 20, 'linewidth', 2, ...
  datOut_444x.fevalCountVals, datOut_444x.fNormVals+epsViz, 'v-', 'markersize', 20, 'linewidth', 2    );
grid on;
ylabel( "||f||" );
xlabel( "feval count" );
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_800.iterCountVals, datOut_800.fNormVals+epsViz, 'x-', 'markersize', 20, 'linewidth', 2, ...
  datOut_900.iterCountVals, datOut_900.fNormVals+epsViz, '+-', 'markersize', 20, 'linewidth', 2, ...
  datOut_940.iterCountVals, datOut_940.fNormVals+epsViz, '*-', 'markersize', 20, 'linewidth', 2, ...
  datOut_940x200.iterCountVals, datOut_940x200.fNormVals+epsViz, 'p-', 'markersize', 20, 'linewidth', 2, ...
  datOut_904.iterCountVals, datOut_904.fNormVals+epsViz, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_904x.iterCountVals, datOut_904x.fNormVals+epsViz, 's-', 'markersize', 20, 'linewidth', 2, ...
  datOut_444.iterCountVals, datOut_444.fNormVals+epsViz, '^-', 'markersize', 20, 'linewidth', 2, ...
  datOut_444x.iterCountVals, datOut_444x.fNormVals+epsViz, 'v-', 'markersize', 20, 'linewidth', 2 );
grid on;
ylabel( "||f||" );
xlabel( "iteration count" );
