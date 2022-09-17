function grootXDatOut = groot_x( funchF, vecX0, algoSetPrm=[], default_solverPrm=[], xPrm=[] )
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
	%
	dateStr = datestr(now,31);
	dateStr(" "==dateStr) = "_";
	dateStr("-"==dateStr) = "";
	dateStr(":"==dateStr) = "";
	xName = [ "A" num2str(algoSetPrm.n) "_X" num2str(sizeX) "_F" num2str(sizeF) "__" dateStr ];
	%
	grootXDatOut.algoSetPrm = algoSetPrm;
	grootXDatOut.default_solverPrm = default_solverPrm;
	grootXDatOut.xPrm = xPrm;
	grootXDatOut.xName = xName;
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
		[ this_vecXBest, this_grootFlag, this_fevalCount, this_datOut ] = this_s.f( funchF, vecX0, this_solverPrm );
		this_elapsedTime = time()-this_startTime();
		assert( isrealarray(this_vecXBest,[sizeX,1]) );
		assert( ischar(this_grootFlag) );
		assert( isscalar(this_grootFlag) );
		assert( isposintscalar(this_fevalCount) );
		this_vecFBest = funchF( this_vecXBest );
		assert( isrealarray(this_vecFBest,[sizeF,1]) );
		this_fBest = norm(this_vecFBest);
		%
		switch ( this_grootFlag )
		case GROOT_FLAG__CNVG
			if ( ~isempty(mygetfield( this_solverPrm, "fTol", [] )) )
				assert( this_fBest <= this_solverPrm.fTol*(1.0+100.0*eps) );
			endif
			if ( ~isempty(mygetfield( this_solverPrm, "fevalLimit", [] )) )
				assert( this_fevalCount <= this_solverPrm.fevalLimit );
			endif
		case { GROOT_FLAG__STOP, GROOT_FLAG__FAIL }
			if ( ~isempty(mygetfield( this_solverPrm, "fTol", [] )) )
			if ( ~isempty(mygetfield( this_solverPrm, "fevalLimit", [] )) )
			if ( this_fBest <= this_solverPrm.fTol*(1.0-100.0*eps) && this_fevalCount <= this_solverPrm.fevalLimit )
				flaggedlog( __FILE__, __LINE__, "*** A SOLVE THAT WAS MAKRED AS NON-CNVG WAS CLEARLY CNVGD! ***" );
			endif
			endif
			endif
		otherwise
			error(["Unsupported value of grootFlag (\"" grootFlag "\")."]);
		endswitch
		%
		if ( xPrm.verbLev >= VERBLEV__PROGRESS )
			msg( __FILE__, __LINE__, sprintf( ...
			  "     Algo%2d,  %15s:  %c,  (%10.3e);   %6d,  (%0.3g)", ...
			  algoIndex, this_strSolverName, ...
			  this_grootFlag, this_fBest, this_fevalCount, this_elapsedTime ) );
		elseif ( xPrm.verbLev >= VERBLEV__INFO )
			printf( "%c", this_grootFlag );
		endif
		%
		grootXDatOut.s(algoIndex).strSolverName = this_strSolverName;
		grootXDatOut.s(algoIndex).fBest = this_fBest;
		grootXDatOut.s(algoIndex).grootFlag = this_grootFlag;
		grootXDatOut.s(algoIndex).fevalCount = this_fevalCount;
		grootXDatOut.s(algoIndex).matInfoA = this_datOut.matInfoA;
		grootXDatOut.s(algoIndex).matInfoB = mygetfield( this_datOut, "matInfoB", [] );
		grootXDatOut.s(algoIndex).elapsedTime = this_elapsedTime;
		%
		clear this_*;
	endfor
	grootXDatOut.elapsedTime = time() - startTime;
return;
endfunction
