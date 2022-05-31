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
sizeX = 50;
fType = 1010;
%fSeed = 81871680;
fSeed = -1;
%fSeed = 30916624;
%fSeed = 21203424;
%fSeed = 42795616;
%fSeed = 12824896;
%fSeed = 73603072;
zcompare__setF;
runFStr = sprintf( "zcompare %d_%dx%d", fType, fSeed, sizeX );
msg( __FILE__, __LINE__, sprintf( "Generated F '%s'.", runFStr ) );
%
if (1)
	msg( __FILE__, __LINE__, "Doing simple one-shot." );
	prm = [];
	prm.verbLev = VERBLEV__COPIOUS;
	prm.valdLev = VALDLEV__VERY_HIGH;
	%prm.iterMax = 2000; zlinsolf100( funchF, vecX0, [], prm );
	%zlinsolf195( funchF, vecX0, [], prm );
	sxsolf100( funchF, vecX0, [] , prm );
	msg( __FILE__, __LINE__, "Goodbye." );
	return;
endif
%
%
%
numRuns = 0;
if (0)
%numRuns++; runList(numRuns).r.runType = 800; runList(numRuns).r.prm.iterMax = 15; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm.iterMax = 900; runList(numRuns).r.prmMemo = "";
%resumeRun = numRuns;
load("dat/vecX0_10_55323424x50_resume.m");
numRuns++; runList(numRuns).r.runType = 800; runList(numRuns).r.prm.iterMax = 100; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm.iterMax = 2000; runList(numRuns).r.prmMemo = "";
resumeRun = 0;
else
%numRuns++; runList(numRuns).r.runType = 550; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 800; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 500; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 3000; runList(numRuns).r.prmMemo = "";
%numRuns++; runList(numRuns).r.runType = 1100; runList(numRuns).r.prm = []; runList(numRuns).r.prm.iterMax = 3000; curveType = "b"; runList(numRuns).r.prmMemo = "curve b";
%numRuns++; runList(numRuns).r.runType = 1195; runList(numRuns).r.prm = []; runList(numRuns).r.prm.verbLev = VERBLEV__DETAILED+10; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 1195; runList(numRuns).r.prm = []; runList(numRuns).r.prmMemo = "";

if (0)
numRuns++; runList(numRuns).r.runType = 1150; runList(numRuns).r.prm.iterMax = 3000; runList(numRuns).r.prmMemo = "";
numRuns++; runList(numRuns).r.runType = 1150; runList(numRuns).r.prm.iterMax = 3000; runList(numRuns).r.prm.epsFD = 1.0E-6; runList(numRuns).r.prmMemo = "eps-6";
numRuns++; runList(numRuns).r.runType = 1150; runList(numRuns).r.prm.iterMax = 3000; runList(numRuns).r.prm.epsFD = 1.0E-5; runList(numRuns).r.prmMemo = "eps-5";
numRuns++; runList(numRuns).r.runType = 1150; runList(numRuns).r.prm.iterMax = 3000; runList(numRuns).r.prm.epsFD = 1.0E-4; runList(numRuns).r.prmMemo = "eps-4";
numRuns++; runList(numRuns).r.runType = 1150; runList(numRuns).r.prm.iterMax = 3000; runList(numRuns).r.prm.epsFD = 1.0E-3; runList(numRuns).r.prmMemo = "eps-3";
numRuns++; runList(numRuns).r.runType = 1150; runList(numRuns).r.prm.iterMax = 3000; runList(numRuns).r.prm.epsFD = 1.0E-2; runList(numRuns).r.prmMemo = "eps-2";
numRuns++; runList(numRuns).r.runType = 1150; runList(numRuns).r.prm.iterMax = 3000; runList(numRuns).r.prm.epsFD = 1.0E-1; runList(numRuns).r.prmMemo = "eps-1";
numRuns++; runList(numRuns).r.runType = 1150; runList(numRuns).r.prm.iterMax = 3000; runList(numRuns).r.prm.epsFD = 1.0E0; runList(numRuns).r.prmMemo = "eps0";
numRuns++; runList(numRuns).r.runType = 1150; runList(numRuns).r.prm.iterMax = 3000; runList(numRuns).r.prm.epsFD = 1.0E+1; runList(numRuns).r.prmMemo = "eps+1";
endif

if (0)
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useQuadUpdate = false;
	runList(numRuns).r.prmMemo = "quad false";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useQuadUpdate = true;
	runList(numRuns).r.prmMemo = "quad true";
endif

if (0)
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useBBall = false;
	runList(numRuns).r.prmMemo = "bball false";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useBBall = true;
	runList(numRuns).r.prmMemo = "bball true";
endif

if (0)
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useQuadUpdate = false;
	runList(numRuns).r.prmMemo = "quad false";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useQuadUpdate = true;
	runList(numRuns).r.prmMemo = "quad true";
endif

if (0)
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useDogLeg = false;
	runList(numRuns).r.prmMemo = "dogLeg false";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useDogLeg = true;
	runList(numRuns).r.prmMemo = "dogLeg true";
endif

if (0)
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useDogLeg = false;
	runList(numRuns).r.prm.curveType = "1";
	runList(numRuns).r.prmMemo = "Lev*1";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useDogLeg = true;
	runList(numRuns).r.prm.curveType = "1";
	runList(numRuns).r.prmMemo = "Powell*1";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useDogLeg = false;
	runList(numRuns).r.prm.curveType = "b";
	runList(numRuns).r.prmMemo = "Lev*BTB";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useDogLeg = true;
	runList(numRuns).r.prm.curveType = "b";
	runList(numRuns).r.prmMemo = "Powell*BTB";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useDogLeg = false;
	runList(numRuns).r.prm.curveType = "m";
	runList(numRuns).r.prmMemo = "Lev*Marq";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useDogLeg = true;
	runList(numRuns).r.prm.curveType = "m";
	runList(numRuns).r.prmMemo = "Powell*Marq";
endif

if (0)
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	[ fooL, fooU ] = lu( matJA );
	runList(numRuns).r.prm.precon_matL = fooL;
	runList(numRuns).r.prm.precon_matU = fooU;
	runList(numRuns).r.prmMemo = "precon:LU";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.precon_matJA = matJA;
	runList(numRuns).r.prmMemo = "precon:JA";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.funchPrecon = @(rhoF,x,f)( matJA \ rhoF );
	runList(numRuns).r.prmMemo = "precon:funch";
endif

if (0)
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useQuadUpdate = false;
	runList(numRuns).r.prm.useBBall = false;
	runList(numRuns).r.prmMemo = "quad- bball-";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useQuadUpdate = true;
	runList(numRuns).r.prm.useBBall = false;
	runList(numRuns).r.prmMemo = "quad+ bball-";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useQuadUpdate = false;
	runList(numRuns).r.prm.useBBall = true;
	runList(numRuns).r.prmMemo = "quad- bball+";
	%
	numRuns++;
	runList(numRuns).r.runType = 1100;
	runList(numRuns).r.prm = [];
	runList(numRuns).r.prm.iterMax = 3000;
	runList(numRuns).r.prm.useQuadUpdate = true;
	runList(numRuns).r.prm.useBBall = true;
	runList(numRuns).r.prmMemo = "quad+ bball+";
endif


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
