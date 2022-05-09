clear;
setVerbLevs;
%setprngstates(0); %500x500 tricky.
%setprngstates(90186240); %For 150x150 triggers levPull in slinsolf200.
%setprngstates(60288928); %200x200 stallz (except 800), triggers path in findZero_904
%matJE = diag(1.0+abs(randn(min([sizeF,sizeX]),1))) + 0.3*randn(sizeF,sizeX); 1.0e-4
%
%setprngstates(67710144); %200x200, 940x200 is bad.
%setprngstates(30223904);% 300x300 fsolve converges; others fail.
%matJE = diag(10.0+abs(randn(min([sizeF,sizeX]),1))) + 0.3*randn(sizeF,sizeX); 1.0e-2
setprngstates(2315152); %300x300: 940x200 good.
%
numFigs = 0;
%
%sizeX = 500; sizeF = 500;
sizeX = 300; sizeF = 300;
%sizeX = 200; sizeF = 200;
%sizeX = 150; sizeF = 150;
%sizeX = 50; sizeF = 50;
%
vecXE = randn(sizeX,1);
matJE = diag(10.0+abs(randn(min([sizeF,sizeX]),1))) + 0.3*randn(sizeF,sizeX);
matA0 = 1.0e-2*randn(sizeF,sizeX);
matA1 = randn(sizeX,sizeX);
matA2 = randn(sizeX,sizeX);
matB0 = 1.0e-2*randn(sizeF,sizeX);
matB1 = randn(sizeX,sizeX);
matB2 = randn(sizeX,sizeX);
matB3 = randn(sizeX,sizeX);
y = @(x)( x - vecXE );
funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
msg( __FILE__, __LINE__, sprintf( "rcond(matJE'*matJE) = %0.3e.", rcond(matJE'*matJE) ) );
%
%vecX0 = zeros(sizeX,1);
vecX0 = vecXE + 1.0e-1*randn(sizeX,1);
vecF0 = funchF(vecX0);
%Df = jacobs( vecX0, funchF );
%msg( __FILE__, __LINE__, sprintf( "rcond(Df'*Df) = %0.3e.", rcond(Df'*Df) ) );



prm = [];
prm.iterMax = 2000;
timeSS = time();
[ vecXF_z100, vecFF_z100, datOut_z100 ] = zlinsolf100( funchF, vecX0, [], prm );
time_z100 = time()-timeSS;
%
%msg( __FILE__, __LINE__, sprintf("Elapsed time = %10.3e.",time()-timeSS) ); msg( __FILE__, __LINE__, "Goodbye!" ); return;


timeSS = time();
[ vecXF_fsolve, vecFF_fsolve, datOut_fsolve, ] = findZero_fsolve( vecX0, funchF );
time_fsolve = time()-timeSS;
msg( __FILE__, __LINE__, "fsolve results..." );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "fsolve", norm(vecFF_fsolve), datOut_fsolve.iterCount, datOut_fsolve.fevalCount, time_fsolve ) );
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









icos_z100 = datOut_z100.iterCountOfSteps;
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
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "z100", norm(vecFF_z100), datOut_z100.stepCount, datOut_z100.fevalCount, time_z100 ) );
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
  datOut_904x.fevalCountVals, datOut_904x.fNormVals+epsViz, 's-', 'markersize', 20, 'linewidth', 2, ...
  datOut_z100.fevalCountVals(icos_z100), datOut_z100.fNormVals(icos_z100)+epsViz, '+-', 'markersize', 20, 'linewidth', 2  );
title( "legend" );
legend( ...
  "800", ...
  "940", ...
  "940x200", ...
  "904", ...
  "904x", ...
  "z100", ...
  "location", "southwest" );
%
numFigs++; figure( numFigs );
semilogy( ...
  datOut_800.fevalCountVals, datOut_800.fNormVals+epsViz, 'x-', 'markersize', 20, 'linewidth', 2, ...
  datOut_940.fevalCountVals, datOut_940.fNormVals+epsViz, '*-', 'markersize', 20, 'linewidth', 2, ...
  datOut_940x200.fevalCountVals, datOut_940x200.fNormVals+epsViz, 'p-', 'markersize', 20, 'linewidth', 2, ...
  datOut_904.fevalCountVals, datOut_904.fNormVals+epsViz, 'o-', 'markersize', 20, 'linewidth', 2, ...
  datOut_904x.fevalCountVals, datOut_904x.fNormVals+epsViz, 's-', 'markersize', 20, 'linewidth', 2, ...
  datOut_z100.fevalCountVals(icos_z100), datOut_z100.fNormVals(icos_z100)+epsViz, '+-', 'markersize', 20, 'linewidth', 2  );
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
  datOut_904x.iterCountVals, datOut_904x.fNormVals+epsViz, 's-', 'markersize', 20, 'linewidth', 2, ...
  datOut_z100.iterCountOfSteps, datOut_z100.fNormVals(icos_z100)+epsViz, '+-', 'markersize', 20, 'linewidth', 2 );
grid on;
ylabel( "||f||" );
xlabel( "step count" );
