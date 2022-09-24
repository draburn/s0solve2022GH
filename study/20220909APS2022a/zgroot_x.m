function grootXDatOut = zgroot_x( zcdo, probIndex=1, prm=[] );
	backup_prngStateDat = getprngstatedat();
	mydefs;
	startTime = time();
	%
	probSetPrm = zcdo.probSetPrm;
	algoSetPrm = zcdo.algoSetPrm;
	default_solverPrm = zcdo.default_solverPrm;
	assert( isposintscalar(probIndex) );
	assert( 1 <= probIndex );
	assert( probIndex <= probSetPrm.numProbs );
	%
	setprngstates( probSetPrm.setSeed, false );
	probSeeds = floor( 1E8*rand(1,probSetPrm.numProbs) );
	setprngstatedat(backup_prngStateDat); % Do this ASAP, in case of abrupt return.
	%
	probType = probSetPrm.probType;
	sizeXIsh = probSetPrm.numUnknowns;
	probSeed = probSeeds(probIndex);
	probGenPrm = [];
	%
	[ funchF, vecX0, probDat ] = genFunchAPS2022_fromType( probType, sizeXIsh, probSeed, probGenPrm );
	assert( is_function_handle(funchF) );
	assert( isrealarray(vecX0,[probDat.sizeX,1]) );
	probStr = sprintf( "Prob %d", probIndex);
	probSizeStr = sprintf( "(%dx%d)", probDat.sizeF, probDat.sizeX );
	msg( __FILE__, __LINE__, sprintf( " %15s  %15s....", probStr, probSizeStr ) );
	%
	xPrm = zcdo.prob(probIndex).xPrm;
	xPrm.verbLev = VERBLEV__PROGRESS;
	default_solverPrm.verbLev = VERBLEV__PROGRESS;
	%
	grootXDatOut = groot_x( funchF, vecX0, algoSetPrm, default_solverPrm, xPrm );
return;
endfunction
