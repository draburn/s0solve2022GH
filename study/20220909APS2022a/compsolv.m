function compsolvDatOut = compsolv( probSetPrm=[], algoSetPrm=[], default_solverPrm=[], compsolvPrm=[] )
	backup_prngStateDat = getprngstatedat();
	mydefs;
	startTime = time();
	if ( isempty(probSetPrm) )
		probSetPrm.numProbs = 10;
		probSetPrm.probType = "test1";
		probSetPrm.numUnknowns = 20;
		probSetPrm.setSeed = 0;
		compsolvPrm.verbLev = mygetfield( compsolvPrm, "verbLev", VERBLEV__INFO );
	endif
	if ( isempty(algoSetPrm) )
		algoSetPrm.n = 3;
		algoSetPrm.s(1).f = @groot_jfnk_basic;
		algoSetPrm.s(1).p.btCoeff = 0.0;
		algoSetPrm.s(2).f = @groot_jfnk_basic;
		algoSetPrm.s(3).f = @groot_fsolve;
	endif
	compsolvPrm.verbLev = mygetfield( compsolvPrm, "verbLev", VERBLEV__WARNING );
	compsolvPrm.valdLev = mygetfield( compsolvPrm, "valdLev", VALDLEV__HIGH );
	%
	msgif( compsolvPrm.verbLev >= VERBLEV__INFO, __FILE__, __LINE__, sprintf( ...
	  "Started zapgroot with numProbs = %d, probType = \"%s\", numUnknowns = %d, setSeed = %d.", ...
	  probSetPrm.numProbs, probSetPrm.probType, probSetPrm.numUnknowns, probSetPrm.setSeed ) );
	%
	%
	setprngstates( probSetPrm.setSeed, false );
	probSeeds = floor( 1E8*rand(1,probSetPrm.numProbs) );
	setprngstatedat(backup_prngStateDat); % Do this ASAP, in case of abrupt return.
	compsolvDatOut.probSeeds = probSeeds;
	%
	assert( isposintscalar(probSetPrm.numProbs) );
	assert( isposintscalar(probSetPrm.numUnknowns) );
	compsolvDatOut.probSetPrm = probSetPrm;
	compsolvDatOut.algoSetPrm = algoSetPrm;
	compsolvDatOut.default_solverPrm = default_solverPrm;
	compsolvDatOut.compsolvPrm = compsolvPrm;
	doCharLog = (  ( compsolvPrm.verbLev >= VERBLEV__INFO ) && ( compsolvPrm.verbLev < VERBLEV__PROGRESS )  );
	for probIndex = 1 : probSetPrm.numProbs
		if ( stopsignalpresent() )
			msg( __FILE__, __LINE__, "Found stopsignal." );
			break;
		endif
		this_probType = probSetPrm.probType;
		this_sizeXIsh = probSetPrm.numUnknowns;
		this_probSeed = probSeeds(probIndex);
		this_probGenPrm = [];
		%
		compsolvDatOut.prob(probIndex).probType = this_probType;
		compsolvDatOut.prob(probIndex).sizeXIsh = this_sizeXIsh;
		compsolvDatOut.prob(probIndex).probSeed = this_probSeed;
		compsolvDatOut.prob(probIndex).probGenPrm = this_probGenPrm;
		[ this_funchF, this_vecX0, this_probDat ] = genFunchAPS2022_fromType( this_probType, this_sizeXIsh, this_probSeed, this_probGenPrm );
		assert( is_function_handle(this_funchF) );
		assert( isrealarray(this_vecX0,[this_probDat.sizeX,1]) );
		this_probStr = sprintf( "Prob %d", probIndex);
		this_probSizeStr = sprintf( "(%dx%d)", this_probDat.sizeF, this_probDat.sizeX );
		if ( doCharLog )
			msgnnl( __FILE__, __LINE__, sprintf( "  %15s  %15s:  ", this_probStr, this_probSizeStr ) );
		elseif ( compsolvPrm.verbLev >= VERBLEV__PROGRESS )
			msg( __FILE__, __LINE__, sprintf( " %15s  %15s....", this_probStr, this_probSizeStr ) );
		endif
		compsolvDatOut.prob(probIndex).sizeX = this_probDat.sizeX;
		compsolvDatOut.prob(probIndex).sizeF = this_probDat.sizeF;
		%
		this_xPrm = [];
		if ( doCharLog )
			this_xPrm.verbLev = VERBLEV__INFO;
		endif
		this_grootXDatOut = groot_x( this_funchF, this_vecX0, algoSetPrm, default_solverPrm, this_xPrm );
		compsolvDatOut.prob(probIndex).grootXDatOut = this_grootXDatOut;
		%
		if ( doCharLog )	
			printf( "\n" );
		endif
		clear this_*;
	endfor
	if ( compsolvPrm.verbLev >= VERBLEV__INFO )
		numAlgos = mygetfield( algoSetPrm, "n", -1 );
		msg( __FILE__, __LINE__, sprintf( "Completed %dx%d algo(solves) in %gs.", probSetPrm.numProbs, numAlgos, time()-startTime )  );
	endif
	compsolvDatOut.elapsedTime = time()-startTime;
	setprngstatedat(backup_prngStateDat); % May be redundant.
return;
endfunction
