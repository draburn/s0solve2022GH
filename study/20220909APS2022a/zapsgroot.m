function zapsgroot = zapgroot( probSetPrm=[], algoSetPrm=[], default_solverPrm=[], zPrm=[] )
	backup_prngStateDat = getprngstatedat();
	mydefs;
	startTime = time();
	if ( isempty(probSetPrm) )
		probSetPrm.numProbs = 3;
		probSetPrm.probType = "test1";
		probSetPrm.numUnknowns = 50;
		probSetPrm.setSeed = 1;
		zPrm.verbLev = mygetfield( zPrm, "verbLev", VERBLEV__INFO );
	endif
	zPrm.verbLev = mygetfield( zPrm, "verbLev", VERBLEV__WARNING );
	zPrm.valdLev = mygetfield( zPrm, "valdLev", VALDLEV__HIGH );
	%
	msgif( zPrm.verbLev >= VERBLEV__INFO, __FILE__, __LINE__, sprintf( ...
	  "Started zapgroot with numProbs = %d, probType = \"%s\", numUnknowns = %d, setSeed = %d.", ...
	  probSetPrm.numProbs, probSetPrm.probType, probSetPrm.numUnknowns, probSetPrm.setSeed ) );
	%
	%
	setprngstates( probSetPrm.setSeed, false );
	probSeeds = floor( 1E8*rand(1,probSetPrm.numProbs) );
	setprngstatedat(backup_prngStateDat); % Do this ASAP, in case of abrupt return.
	%
	xPrm = [];
	%
	assert( isposintscalar(probSetPrm.numProbs) );
	assert( isposintscalar(probSetPrm.numUnknowns) );
	zDatOut.probSetPrm = probSetPrm;
	zDatOut.algoSetPrm = algoSetPrm;
	zDatOut.default_solverPrm = default_solverPrm;
	zDatOut.zPrm = zPrm;
	zDatOut.xPrm = xPrm;
	doCharLog = (  ( zPrm.verbLev >= VERBLEV__INFO ) && ( zPrm.verbLev < VERBLEV__PROGRESS )  );
	for probIndex = 1 : probSetPrm.numProbs
		if ( stopsignalpresent() )
			msg( __FILE__, __LINE__, "Found stopsignal." );
			break;
		endif
		this_probSeed = probSeeds(probIndex);
		this_probGenPrm = [];
		[ this_funchF, this_vecX0, this_probDat ] = genFunchAPS2022_fromType( probSetPrm.probType, probSetPrm.numUnknowns, this_probSeed, this_probGenPrm );
		assert( is_function_handle(this_funchF) );
		assert( isrealarray(this_vecX0,[this_probDat.sizeX,1]) );
		this_probStr = sprintf( "Prob %d", probIndex);
		this_probSizeStr = sprintf( "(%dx%d)", this_probDat.sizeF, this_probDat.sizeX );
		if ( doCharLog )
			msgnnl( __FILE__, __LINE__, sprintf( "  %15s  %s:  ", this_probStr, this_probSizeStr ) );
		elseif ( zPrm.verbLev >= VERBLEV__PROGRESS )
			msg( __FILE__, __LINE__, sprintf( " %15s  %s....", this_probStr, this_probSizeStr ) );
		endif
		%
		zDatOut.prob(probIndex).probSeed = this_probSeed;
		zDatOut.prob(probIndex).probGenPrm = this_probGenPrm;
		zDatOut.prob(probIndex).sizeX = this_probDat.sizeX;
		zDatOut.prob(probIndex).sizeF = this_probDat.sizeF;
		%
		this_xPrm = [];
		if ( doCharLog )
			this_xPrm.verbLev = VERBLEV__INFO;
		endif
		this_grootxDatOut = groot_x( this_funchF, this_vecX0, algoSetPrm, default_solverPrm, this_xPrm );
		%
		if ( doCharLog )	
			printf( "\n" );
		endif
		clear this_*;
	endfor
	zDatOut.probSeeds = probSeeds;
	setprngstatedat(backup_prngStateDat); % May be redundant.
return;
endfunction
