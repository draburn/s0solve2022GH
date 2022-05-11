clear;
setVerbLevs;
% This is about as hard as I can make this model and still reliably get convergence.
%
setprngstates();
%
numFigs = 0;
%
sizeX = 1000; sizeF = 1000;
%sizeX = 500; sizeF = 500;
%sizeX = 300; sizeF = 300;
%sizeX = 200; sizeF = 200;
%sizeX = 150; sizeF = 150;
%sizeX = 50; sizeF = 50;
%
myrand = @(s1,s2)( 2.0*rand(s1,s2) - 1.0 );
%
vecXE = myrand(sizeX,1);
matJE = eye(sizeF,sizeX) + 1.0e-2*myrand(sizeF,sizeX);
matA0 = 0*1.0e-4*myrand(sizeF,sizeX);
matA1 = myrand(sizeX,sizeX);
matA2 = myrand(sizeX,sizeX);
matB0 = 0*1.0e-5*myrand(sizeF,sizeX);
matB1 = myrand(sizeX,sizeX);
matB2 = myrand(sizeX,sizeX);
matB3 = myrand(sizeX,sizeX);
y = @(x)( x - vecXE );
funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
msg( __FILE__, __LINE__, sprintf( "rcond(matJE'*matJE) = %0.3e.", rcond(matJE'*matJE) ) );
%
vecX0 = myrand(sizeX,1);
vecF0 = funchF(vecX0);
%Df = jacobs( vecX0, funchF );
%msg( __FILE__, __LINE__, sprintf( "rcond(Df'*Df) = %0.3e.", rcond(Df'*Df) ) );


if (0)
prm = [];
prm.iterMax = 2000;
timeSS = time();
[ vecXF_z100, vecFF_z100, datOut_z100 ] = zlinsolf100( funchF, vecX0, [], prm );
time_z100 = time()-timeSS;
endif
%
%msg( __FILE__, __LINE__, sprintf("Elapsed time = %10.3e.",time()-timeSS) ); msg( __FILE__, __LINE__, "Goodbye!" ); return;


timeSS = time();
[ vecXF_fsolve, vecFF_fsolve, datOut_fsolve, ] = findZero_fsolve( vecX0, funchF );
time_fsolve = time()-timeSS;
msg( __FILE__, __LINE__, "fsolve results..." );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "fsolve", norm(vecFF_fsolve), datOut_fsolve.iterCount, datOut_fsolve.fevalCount, time_fsolve ) );
if ( norm(vecFF_fsolve) > sqrt(eps) )
	msg( __FILE__, __LINE__, "fsolve() failed to converge. Goodbye!" ); return;
endif
%msg( __FILE__, __LINE__, "Goodbye!" ); return;

%msg( __FILE__, __LINE__, "Goodbye!" ); return;



prm = [];
prm.iterMax = 20;
prm.slinsolfver = 100;
timeSS = time();
[ vecXF_940, vecFF_940, datOut_940 ] = findZero_940( vecX0, funchF, prm );
time_940 = time()-timeSS;
%
%
%
prm = [];
prm.step_prm.usePostLinsolfPhiPatch = false;
prm.iterMax = 20;
timeSS = time();
[ vecXF_904, vecFF_904, datOut_904 ] = findZero_904( vecX0, funchF, prm );
time_904 = time()-timeSS;
%
%
prm = [];
prm.step_prm.usePostLinsolfPhiPatch = true;
prm.iterMax = 20;
timeSS = time();
[ vecXF_904x, vecFF_904x, datOut_904x ] = findZero_904( vecX0, funchF, prm );
time_904x = time()-timeSS;
%
%
prm = [];
prm.iterMax = 20;
timeSS = time();
[ vecXF_800, vecFF_800, datOut_800 ] = findZero_800( vecX0, funchF, prm );
time_800 = time()-timeSS;


prm = [];
prm.iterMax = 20;
prm.slinsolfver = 200;
timeSS = time();
[ vecXF_940x200, vecFF_940x200, datOut_940x200 ] = findZero_940( vecX0, funchF, prm );
time_940x200 = time()-timeSS;
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "940x200", norm(vecFF_940x200), datOut_940x200.iterCount, datOut_940x200.fevalCount, time_940x200 ) );


%if (stopsignalpresent())
%	msg( __FILE__, __LINE__, "Received stop signal." );
%	return;
%endif
%
%









%icos_z100 = datOut_z100.iterCountOfSteps;
% Push to triple-digit line number.
msg( __FILE__, __LINE__, sprintf( "norm(vecF0) = %g.", norm(vecF0) ) );
msg( __FILE__, __LINE__, "Solver results..." );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11s;  %11s;  %11s;  %11s.", "solver name", "||vecFF||", "stepCount", "fevalCount", "time(s)" ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "fsolve", norm(vecFF_fsolve), datOut_fsolve.iterCount, datOut_fsolve.fevalCount, time_fsolve ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "800", norm(vecFF_800), datOut_800.iterCount, datOut_800.fevalCount, time_800 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "940", norm(vecFF_940), datOut_940.iterCount, datOut_940.fevalCount, time_940 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "940x200", norm(vecFF_940x200), datOut_940x200.iterCount, datOut_940x200.fevalCount, time_940x200 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "904", norm(vecFF_904), datOut_904.iterCount, datOut_904.fevalCount, time_904 ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "904x", norm(vecFF_904x), datOut_904x.iterCount, datOut_904x.fevalCount, time_904x ) );
%msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "z100", norm(vecFF_z100), datOut_z100.iterCount, datOut_z100.fevalCount, time_z100 ) );
%msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "z100", norm(vecFF_z100), datOut_z100.stepCount, datOut_z100.fevalCount, time_z100 ) );
%
%
%
epsViz = 1.0e-18;
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_800.fevalCountVals, datOut_800.fNormVals+epsViz, 'x-', 'markersize', 20, 'linewidth', 2, ...
  datOut_940.fevalCountVals, datOut_940.fNormVals+epsViz, '*-', 'markersize', 20, 'linewidth', 2, ...
  datOut_940x200.fevalCountVals, datOut_940x200.fNormVals+epsViz, 'p-', 'markersize', 20, 'linewidth', 2, ...
  datOut_904.fevalCountVals, datOut_904.fNormVals+epsViz, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_904x.fevalCountVals, datOut_904x.fNormVals+epsViz, 's-', 'markersize', 20, 'linewidth', 2  );
title( "legend" );
legend( ...
  "800", ...
  "940", ...
  "940x200", ...
  "904", ...
  "904x", ...
  "location", "northeast" );
grid on;
msg( __FILE__, __LINE__, "Goodbye." ); return;
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_800.fevalCountVals, datOut_800.fNormVals+epsViz, 'x-', 'markersize', 20, 'linewidth', 2, ...
  datOut_940.fevalCountVals, datOut_940.fNormVals+epsViz, '*-', 'markersize', 20, 'linewidth', 2, ...
  datOut_940x200.fevalCountVals, datOut_940x200.fNormVals+epsViz, 'p-', 'markersize', 20, 'linewidth', 2, ...
  datOut_904.fevalCountVals, datOut_904.fNormVals+epsViz, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_904x.fevalCountVals, datOut_904x.fNormVals+epsViz, 's-', 'markersize', 20, 'linewidth', 2  );
grid on;
ylabel( "||f||" );
xlabel( "feval count" );
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_800.iterCountVals, datOut_800.fNormVals+epsViz, 'x-', 'markersize', 20, 'linewidth', 2, ...
  datOut_940.iterCountVals, datOut_940.fNormVals+epsViz, '*-', 'markersize', 20, 'linewidth', 2, ...
  datOut_940x200.iterCountVals, datOut_940x200.fNormVals+epsViz, 'p-', 'markersize', 20, 'linewidth', 2, ...
  datOut_904.iterCountVals, datOut_904.fNormVals+epsViz, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_904x.iterCountVals, datOut_904x.fNormVals+epsViz, 's-', 'markersize', 20, 'linewidth', 2 );
grid on;
ylabel( "||f||" );
xlabel( "step count" );
