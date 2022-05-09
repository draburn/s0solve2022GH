clear;
setVerbLevs;
numFigs = 0;
mainStartTime = time();
mainStartDatestr = datestr(now,31);
msg( __FILE__, __LINE__, sprintf( "Starting run suite '%s'.", mainStartDatestr ) );
%
%
sizeX = 10;
fType = 10;
fSeed = 0;
zcompare__setF;
%
%
%
numRuns = 0;
numRuns++; runList(numRuns).r.runType = 50; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 550; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "(WIP)";
numRuns++; runList(numRuns).r.runType = -1; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "(dne)";
%
%
assert( 0 < numRuns );
for runIndex = 1 : numRuns
	r = runList(runIndex).r;
	if (isempty(r.prmMemo))
		r.runName = sprintf( "%d", r.runType );
	else
		r.runName = sprintf( "%d %s", r.runType, r.prmMemo );
	endif
	msg( __FILE__, __LINE__, sprintf( "Starting run %d/%d: '%s'...", runIndex, numRuns, r.runName ) );
	zcompare__doRun;
	markerTypes = "+o*xsd^v<>ph";
	r.mlStyle = [ markerTypes(mod(runIndex,numRuns)+1), "-" ];
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
		semilogy( r.fevalCountOfStep, r.fBestNormOfStep+epsViz, r.mlStyle, 'markersize', 20, 'linewidth', 2 );
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
		semilogy( r.fevalCountOfStep, r.fBestNormOfStep+epsViz, r.mlStyle, 'markersize', 20, 'linewidth', 2 );
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
		semilogy( (0:r.stepCount), r.fBestNormOfStep+epsViz, r.mlStyle, 'markersize', 20, 'linewidth', 2 );
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
