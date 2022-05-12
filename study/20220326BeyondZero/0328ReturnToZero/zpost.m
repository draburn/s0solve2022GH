plotStartTime = time();
vecF0 = funchF(vecX0);
%
msg( __FILE__, __LINE__, sprintf( "Run suite '%s' with F '%s' completed in %0.3es.", mainStartDatestr, runFStr, mainCalcElapsedTime ) );
%
msg( __FILE__, __LINE__, sprintf( "norm(vecF0) = %g.", norm(vecF0) ) );
msg( __FILE__, __LINE__, "Solver results..." );
msg( __FILE__, __LINE__, sprintf( " %15s:  %11s;  %11s;  %11s;  %11s.", "solver name", "||vecFF||", "stepCount", "fevalCount", "time(s)" ) );
for runIndex = 1 : numRuns
	r = runList(runIndex).r;
	if (r.isValid)
		msg( __FILE__, __LINE__, sprintf( " %15s:  %11.3e;  %11d;  %11d;  %11.3e.", ...
		  r.runName, norm(r.vecFF), r.stepCount, r.fevalCount, r.elapsedTime ) );
	else
		msg( __FILE__, __LINE__, sprintf( " %15s:  %11.3e;  %11d;  %11d;  %11.3e.", r.runName, -1.0, -1, -1, -1.0 ) );
	endif
endfor
%
%
msg( __FILE__, __LINE__, "Generating plots..." );
%
%
numFigs++; figure(numFigs);
epsViz = 1.0e-18;
leg = {};
m = 0;
for n=1:numRuns
	r = runList(n).r;
	if (r.isValid)
		%loglog( r.fevalCountOfStep, r.fBestNormOfStep+epsViz, r.mlStyle, 'markersize', r.mSize, 'linewidth', 2 );
		semilogy( r.fevalCountOfStep, r.fBestNormOfStep+epsViz, r.mlStyle, 'markersize', r.mSize, 'linewidth', 2 );
		hold on;
		m++;
		leg(m) = sprintf( "%d: %s (%s)", n, r.runName, r.runTypeDescrip );
	endif
endfor
hold off;
grid on;
set( legend( leg ), 'Interpreter', 'none' );
set( xlabel( "feval count" ), 'Interpreter', 'none' );
set( ylabel( "||F best||" ), 'Interpreter', 'none' );
set( title([ mainStartDatestr " " runFStr " CNVG V FEVAL" ]), 'Interpreter', 'none' );
%
%
if (1)
numFigs++; figure(numFigs);
epsViz = 1.0e-18;
leg = {};
for n=1:numRuns
	r = runList(n).r;
	if (r.isValid)
		%loglog( r.fevalCountOfStep, r.fBestNormOfStep+epsViz, r.mlStyle, 'markersize', r.mSize, 'linewidth', 2 );
		semilogy( r.fevalCountOfStep, r.fBestNormOfStep+epsViz, r.mlStyle, 'markersize', r.mSize, 'linewidth', 2 );
		hold on;
	endif
endfor
hold off;
grid on;
set( xlabel( "feval count" ), 'Interpreter', 'none' );
set( ylabel( "||F best||" ), 'Interpreter', 'none' );
set( title([ mainStartDatestr " " runFStr " CNVG V FEVAL" ]), 'Interpreter', 'none' );
endif
%
%
if (1)
numFigs++; figure(numFigs);
epsViz = 1.0e-18;
for n=1:numRuns
	r = runList(n).r;
	if (r.isValid)
		%loglog( r.stepCountOfStep, r.fBestNormOfStep+epsViz, r.mlStyle, 'markersize', r.mSize, 'linewidth', 2 );
		semilogy( r.stepCountOfStep, r.fBestNormOfStep+epsViz, r.mlStyle, 'markersize', r.mSize, 'linewidth', 2 );
		hold on;
	endif
endfor
hold off;
grid on;
%set( legend( leg ), 'Interpreter', 'none' );
set( xlabel( "step count" ), 'Interpreter', 'none' );
set( ylabel( "||F best||" ), 'Interpreter', 'none' );
set( title([ mainStartDatestr " " runFStr " CNVG V STEP" ]), 'Interpreter', 'none' );
endif
%
%
mainCPPElapsedTime = time()-mainStartTime;

msg( __FILE__, __LINE__, sprintf( "Run suite '%s' with F '%s' completed in %0.3es; post-proc took %0.3es.", ...
  mainStartDatestr, runFStr, mainCalcElapsedTime, time()-plotStartTime ) );
msg( __FILE__, __LINE__, "Goodbye!" );
printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n" );
