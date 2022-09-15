function xDatOut = groot_x( funchF, vecX0, algoSetPrm=[], default_solverPrm=[], xPrm=[] )
	mydefs;
	startTime = time();
	if ( isempty(algoSetPrm) )
		algoSetPrm.n = 3;
		algoSetPrm.s(1).f = @groot_fsolve;
		algoSetPrm.s(2).f = @groot_jfnk_basic;
		algoSetPrm.s(3).f = @groot_jfnk_basic;
		algoSetPrm.s(3).p.btCoeff = 0.0;
		xPrm.verbLev = mygetfield( xPrm, "verbLev", VERBLEV__PROGRESS );
	endif
	xPrm.verbLev = mygetfield( xPrm, "verbLev", VERBLEV__WARNING );
	xPrm.valdLev = mygetfield( xPrm, "valdLev", VALDLEV__HIGH );
	%
	assert( is_function_handle(funchF) );
	sizeX = size( vecX0, 1 );
	assert( isrealarray(vecX0,[sizeX,1]) );
	vecF0 = funchF( vecX0 );
	sizeF = size( vecF0, 1 );
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	assert( isposintscalar(algoSetPrm.n) );
	assert( isvector(algoSetPrm.s) );
	assert( max(size(algoSetPrm.s)) == algoSetPrm.n );
	for algoIndex = 1 : algoSetPrm.n
		if ( stopsignalpresent() )
			msg( __FILE__, __LINE__, "Found stopsignal." );
			break;
		endif
		this_s = algoSetPrm.s(algoIndex);
		%
		assert( is_function_handle( this_s.f ) );
		this_fInfo = functions( this_s.f );
		this_strSolverName = trimFileName(this_fInfo.function)(7:end);
		if ( ~isempty(mygetfield( this_s, "p", [] )) )
			this_strSolverName = [ "~" this_strSolverName ];
		endif
		%
		this_solverPrm = mygetfield( this_s, "p", [] );
		this_solverPrm = overwritefields( default_solverPrm, this_solverPrm );
		%
		this_startTime = time();
		[ this_vecXBest, this_strGrootFlag, this_fevalCount, this_datOut ] = this_s.f( funchF, vecX0, this_solverPrm );
		this_elapsedTime = time()-this_startTime();
		assert( isrealarray(this_vecXBest,[sizeX,1]) );
		assert( ischar(this_strGrootFlag) );
		assert( isvector(this_strGrootFlag) );
		assert( isposintscalar(this_fevalCount) );
		%
		this_vecFBest = funchF( this_vecXBest );
		assert( isrealarray(this_vecFBest,[sizeF,1]) );
		this_fBest = norm(this_vecFBest);
		%
		if ( xPrm.verbLev >= VERBLEV__PROGRESS )
			msg( __FILE__, __LINE__, sprintf( ...
			  "     Algo%2d,  %15s:  %c,  (%10.3e);   %6d,  (%0.3g)", ...
			  algoIndex, this_strSolverName, ...
			  this_strGrootFlag, this_fBest, this_fevalCount, this_elapsedTime ) );
		elseif ( xPrm.verbLev >= VERBLEV__INFO )
			printf( "%c", this_strGrootFlag );
		endif
		%
		xDatOut.s(algoIndex).fBest = this_fBest;
		xDatOut.s(algoIndex).strGrootFlag = this_strGrootFlag;
		xDatOut.s(algoIndex).fevalCount = this_fevalCount;
		xDatOut.s(algoIndex).matInfoA = this_datOut.matInfoA;
		xDatOut.s(algoIndex).matInfoB = mygetfield( this_datOut, "matInfoB", [] );
		xDatOut.s(algoIndex).elapsedTime = this_elapsedTime;
		%
		clear this_*;
	endfor
	xDatOut.elapsedTime = time() - startTime;
return;
endfunction
