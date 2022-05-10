clear;
setVerbLevs;
numFigs = 0;
mainStartTime = time();
mainStartDatestr = datestr(now,31);
msg( __FILE__, __LINE__, sprintf( "Starting run suite '%s'.", mainStartDatestr ) );
if ( stopsignalpresent() )
	error( "Stop signal is already present." );
endif
%
%
sizeX = 50;
fType = 5;
fSeed = 0;
zcompare__setF;
numRuns = 0;
%
%
%
%numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "(WIP)";
numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm.iterMax = 1000; runList(numRuns).r.prmMemo = "it1000";
numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm.iterMax = 1000; runList(numRuns).r.prmMemo = "it1000";
numRuns++; runList(numRuns).r.runType = 50; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 444; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 550; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 800; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 904; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 904; runList(numRuns).r.prm.slinsolfver = 100; runList(numRuns).r.prmMemo = "sl100";
%numRuns++; runList(numRuns).r.runType = 904; runList(numRuns).r.prm.slinsolfver = 200; runList(numRuns).r.prmMemo = "sl200";
numRuns++; runList(numRuns).r.runType = 940; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 940; runList(numRuns).r.prm.step_prm.usePostLinsolfPhiPatch = false; runList(numRuns).r.prmMemo = "uplpp false";
%numRuns++; runList(numRuns).r.runType = 940; runList(numRuns).r.prm.step_prm.usePostLinsolfPhiPatch = true; runList(numRuns).r.prmMemo = "uplpp true";
%numRuns++; runList(numRuns).r.runType = -1; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "(dne)";
%
%
% Hacks
numRuns = 1;
useResume = true;
haveResumed = false;
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
	r = runList(runIndex).r;
	zcompare__doRun;
	if (useResume)
		if (~haveResumed)
			vecX0 = r.vecXF;
			resumeFevalCount = r.fevalCount;
			resumeStepCount = r.stepCount;
			haveResumed = true;
		else
			r.fevalCountOfStep += resumeFevalCount;
			r.stepCountOfStep += resumeStepCount;
		endif
	endif
	runList(runIndex).r = r;
endfor
vecF0 = funchF(vecX0);
%
mainElapsedTime = time()-runStartTime;
msg( __FILE__, __LINE__, sprintf( "Run suite '%s' with F '%s' completed in %0.3es.", mainStartDatestr, runFStr, mainElapsedTime ) );
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
		leg(m) = sprintf( "%s (%s)", r.runName, r.runTypeDescrip );
	endif
endfor
hold off;
legend( leg );
title( "LEGEND" );
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
msg( __FILE__, __LINE__, "Goodbye!" );
return;
