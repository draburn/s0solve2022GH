clear;
setVerbLevs;
numFigs = 0;
mainStartTime = time();
mainStartDatestr = datestr(now,31);
printf( "\n\nVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n" );
msg( __FILE__, __LINE__, sprintf( "Starting run suite '%s'.", mainStartDatestr ) );
if ( stopsignalpresent() )
	error( "Stop signal is already present." );
endif
%
%
sizeX = 50;
fType = 10;
fSeed = 55323424;
zcompare__setF;
runFStr = sprintf( "zcompare %d_%dx%d", fType, fSeed, sizeX );
msg( __FILE__, __LINE__, sprintf( "Generated F '%s'.", runFStr ) );
%
if (0)
	numRuns = 0;
	runIndex = 0;
	r.runType = 1100; r.prm.iterMax = 1000; r.prmMemo = "it1000";
	zcompare__doRun;
	msg( __FILE__, __LINE__, "Goodbye!" ); return;
endif
%
%
%
numRuns = 0;
%numRuns++; runList(numRuns).r.runType = 550; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 800; runList(numRuns).r.prm.iterMax = 2000; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm.iterMax = 5000; runList(numRuns).r.prmMemo = "";
resumeRun = numRuns;
%
%numRuns++; runList(numRuns).r.runType = 550; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 800; runList(numRuns).r.prm.iterMax = 100; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm.iterMax = 1000; runList(numRuns).r.prmMemo = "";
%
%numRuns++; runList(numRuns).r.runType = 800; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 904; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 904; runList(numRuns).r.prm.slinsolfver = 100; runList(numRuns).r.prmMemo = "sl100";
%numRuns++; runList(numRuns).r.runType = 904; runList(numRuns).r.prm.slinsolfver = 200; runList(numRuns).r.prmMemo = "sl200";
%numRuns++; runList(numRuns).r.runType = 940; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 940; runList(numRuns).r.prm.step_prm.usePostLinsolfPhiPatch = false; runList(numRuns).r.prmMemo = "uplpp false";
%numRuns++; runList(numRuns).r.runType = 940; runList(numRuns).r.prm.step_prm.usePostLinsolfPhiPatch = true; runList(numRuns).r.prmMemo = "uplpp true";
%numRuns++; runList(numRuns).r.runType = -1; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "(dne)";
%
%
%
msg( __FILE__, __LINE__, "TODO: Run with fsolve for comparison." );
%
%
assert( 0 < numRuns );
for runIndex = 1 : numRuns
	if ( stopsignalpresent() )
		msg( __FILE__, __LINE__, "Received stop signal." );
		return;
	endif
	msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
	r = runList(runIndex).r;
	zcompare__doRun;
	if ( 0 < resumeRun )
	if ( runIndex == resumeRun )
		vecX0 = r.vecXF;
		resumeFevalCount = r.fevalCount;
		resumeStepCount = r.stepCount;
	elseif ( runIndex > resumeRun )
		r.fevalCountOfStep += resumeFevalCount;
		r.stepCountOfStep += resumeStepCount;
	endif
	endif
	runList(runIndex).r = r;
	msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
endfor
vecF0 = funchF(vecX0);
%
mainCalcElapsedTime = time()-mainStartTime;
msg( __FILE__, __LINE__, sprintf( "Run suite '%s' with F '%s' completed in %0.3es.", mainStartDatestr, runFStr, mainCalcElapsedTime ) );
%
msg( __FILE__, __LINE__, sprintf( "norm(vecF0) = %g.", norm(vecF0) ) );
msg( __FILE__, __LINE__, "Solver results..." );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11s;  %11s;  %11s;  %11s.", "solver name", "||vecFF||", "stepCount", "fevalCount", "time(s)" ) );
for runIndex = 1 : numRuns
	r = runList(runIndex).r;
	if (r.isValid)
		msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", ...
		  r.runName, norm(r.vecFF), r.stepCount, r.fevalCount, r.elapsedTime ) );
	else
		msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", r.runName, -1.0, -1, -1, -1.0 ) );
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
legObj = legend( leg );
xlabObj = xlabel( "feval count" );
ylabObj = ylabel( "||F best||" );
%title( "||F best|| vs feval count" );
%title( "LEGEND" );
titleObj = title( [ mainStartDatestr " " runFStr " cnvg" ] );
set( legObj, 'Interpreter', 'none' );
set( xlabObj, 'Interpreter', 'none' );
set( ylabObj, 'Interpreter', 'none' );
set( titleObj, 'Interpreter', 'none' );

mainCPPElapsedTime = time()-mainStartTime;
msg( __FILE__, __LINE__, sprintf( "With plots, run suite '%s' with F '%s' completed in %0.3es.", mainStartDatestr, runFStr, time()-mainStartTime ) );
msg( __FILE__, __LINE__, "Goodbye!" );
printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n" );
return;

%
numFigs++; figure(numFigs);
epsViz = 1.0e-18;
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
xlabel( "feval count" );
ylabel( "||F best||" );
title( "||F best|| vs feval count" );
%
%
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
xlabel( "step count" );
ylabel( "||F best||" );
title( "||F best|| vs step count" );
%

mainCPPElapsedTime = time()-mainStartTime;
msg( __FILE__, __LINE__, sprintf( "With plots, run suite '%s' with F '%s' completed in %0.3es.", mainStartDatestr, runFStr, time()-mainStartTime ) );
msg( __FILE__, __LINE__, "Goodbye!" );
return;
