clear;
%%%setVerbLevs;
mydefs;
mainStartTime = time();
mainStartDatestr = datestr(now,31);
printf( "\n\n" );
msg( __FILE__, __LINE__, "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV" );
msg( __FILE__, __LINE__, sprintf( "Starting run suite '%s'.", mainStartDatestr ) );
if ( stopsignalpresent() )
	error( "Stop signal is already present." );
endif
%
%
fType = 1010
fSeed = -1
sizeX = 100
zcompare__setF;
runFStr = sprintf( "zcompare %d_%dx%d", fType, fSeed, sizeX );
msg( __FILE__, __LINE__, sprintf( "Generated F '%s'.", runFStr ) );
%
if (0)
	msg( __FILE__, __LINE__, "Doing simple one-shot." );
	prm = [];
	%prm.verbLev = VERBLEV__DETAILED;
	%prm.verbLev = VERBLEV__UNLIMITED;
	%prm.valdLev = VALDLEV__VERY_HIGH;
	prm.valdLev = VALDLEV__UNLIMITED;
	prm.iterMax = 5000;
	prm.fevalMax = 5000;
	%zlinsolf195( funchF, vecX0, [], prm );
	%sxsolf100( funchF, vecX0, [] , prm );
	%prm.useDogLeg = true;
	%[ vecXF, vecFF, retCode, fevalCount, stepsCount, datOut ] = sxsolf181perCompareAP( funchF, vecX0, [] , prm );
	findZero_800sssl( vecX0, funchF, prm );
	%
	mainCalcElapsedTime = time()-mainStartTime;
	msg( __FILE__, __LINE__, sprintf( "One-shot '%s' with F '%s' completed in %0.3es.", mainStartDatestr, runFStr, mainCalcElapsedTime ) );
	msg( __FILE__, __LINE__, "Goodbye." );
	return;
endif
%
%
%
numRuns = 0;
resumeRun = -1;
%numRuns++; runList(numRuns).r.runType = 800; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 2; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 2181; runList(numRuns).r.prm = []; runList(numRuns).r.prm.stepsMax = 2; runList(numRuns).r.prmMemo = "";

resumeRun = numRuns;

numRuns++; runList(numRuns).r.runType = 800; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 802; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 85551; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 815353; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 810; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 820; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 1000; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 2100; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 1000; runList(numRuns).r.prm.fevalMax = 1000; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 2181; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 1000; runList(numRuns).r.prm.fevalMax = 1000; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 2182; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 1000; runList(numRuns).r.prm.fevalMax = 1000; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 2183; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 1000; runList(numRuns).r.prm.fevalMax = 1000; runList(numRuns).r.prmMemo = "";
%umRuns++; runList(numRuns).r.runType = 21801; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 1000; runList(numRuns).r.prm.fevalMax = 1000; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 21802; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 1000; runList(numRuns).r.prm.fevalMax = 1000; runList(numRuns).r.prmMemo = "";

%numRuns++; runList(numRuns).r.runType = 2124; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 5000; runList(numRuns).r.prm.fevalMax = 5000; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 2114; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 1000; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 1195; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 1000; runList(numRuns).r.prmMemo = "";
%resumeRun = numRuns;
%
%
assert( 0 < numRuns );
for runIndex = 1 : numRuns
	msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
	r = runList(runIndex).r;
	zcompare__doRun;
	if ( stopsignalpresent() )
		msg( __FILE__, __LINE__, "Received stop signal." );
		msg( __FILE__, __LINE__, sprintf( "Run suite '%s' with F '%s' stopped after %0.3es.", ...
		  mainStartDatestr, runFStr, time()-mainStartTime ) );
		msg( __FILE__, __LINE__, "Goodbye!" );
		msg( __FILE__, __LINE__, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" );
		printf("\n\n" );
		return;
	endif
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
doExtras = false;
zpost;
doExtras = true;
