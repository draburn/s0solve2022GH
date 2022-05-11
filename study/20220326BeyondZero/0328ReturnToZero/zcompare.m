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
if (1)
%numRuns++; runList(numRuns).r.runType = 800; runList(numRuns).r.prm.iterMax = 15; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm.iterMax = 900; runList(numRuns).r.prmMemo = "";
%resumeRun = numRuns;
load vecX0.m;
numRuns++; runList(numRuns).r.runType = 800; runList(numRuns).r.prm.iterMax = 100; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm.iterMax = 2000; runList(numRuns).r.prmMemo = "";
resumeRun = 0;
else
%numRuns++; runList(numRuns).r.runType = 550; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm.iterMax = 5000; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 800; runList(numRuns).r.prm.iterMax = 2000; runList(numRuns).r.prmMemo = "";
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
endif
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
mainCalcElapsedTime = time()-mainStartTime;
zpost;
