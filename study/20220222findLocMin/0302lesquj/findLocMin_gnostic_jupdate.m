% Function...
%  A 'gnostic for studying how various forms of Jacobian updating impact convergence.

function [ vecX, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchF, prm=[] )
	commondefs;
	findLocMin_gnostic_jupdate__defs;
	time0 = time();
	fevalCount = 0;
	jevalCount = 0;
	%
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	valdLev = mygetfield( prm, "valdLev", VALDLEV__UNLIMITED );
	if ( valdLev >= VALLEV__LOW )
		assert( isrealscalar(verbLev) );
		assert( isrealscalar(valdLev) );
	endif
	%
	sizeX = size(vecX0,1);
	if ( valdLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecX0,[sizeX,1]) );
	endif
	%
	vecF0 = funchF( vecX0 ); fevalCount++;
	sizeF = size(vecF0,1);
	if ( valdLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecF0,[sizeF,1]) );
	endif
	%
	matJ0 = jacobs( vecX0, funchF ); jevalCount++;
	if ( valdLev >= VALLEV__MEDIUM )
		assert( isrealarray(matJ0,[sizeF,sizeX]) );
	endif
	%
	%
	%
	iterMax = mygetfield( prm, "iterMax", 300 );
	if ( valdLev >= VALLEV__LOW )
		assert( isrealscalar(iterMax) );
	endif
	stepType = mygetfield( prm, "stepType", STEP_TYPE__BLIND_NEWTON );
	if ( valdLev >= VALLEV__LOW )
		assert( isrealscalar(stepType) );
	endif
	%
	%jupdateType = mygetfield( prm, "jupdateType", JUPDATE_TYPE__NONE ); % THIS WORKS BETTER???
	%jupdateType = mygetfield( prm, "jupdateType", JUPDATE_TYPE__BROYDEN );
	jupdateType = mygetfield( prm, "jupdateType", JUPDATE_TYPE__RECALC );
	if ( valdLev >= VALLEV__LOW )
		assert( isrealscalar(jupdateType) );
	endif
	%
	%
	%	
	if ( nargout >= 2 )
		datOut = [];
	endif
	%
	%
	%
	vecX = vecX0;
	vecF = vecF0;
	matJ = matJ0;
	iterCount = 0;
	while (1)
		if ( verbLev >= VERBLEV__PROGRESS )
		if ( 0 == iterCount )
			msg( __FILE__, __LINE__, sprintf( "  %10.3e, %4d;  %5d,  %3d;  %10.3e, %10.3e.", ...
			  time()-time0, iterCount, ...
			  fevalCount, jevalCount, ...
			  sumsq(vecF)/2.0, -1.0 ) );
		else
			msg( __FILE__, __LINE__, sprintf( "  %10.3e, %4d;  %5d,  %3d;  %10.3e, %10.3e.", ...
			  time()-time0, iterCount, ...
			  fevalCount, jevalCount, ...
			  sumsq(vecF)/2.0, (sumsq(vecF_prev)-sumsq(vecF))/2.0 ) );
		endif
		endif
		%
		if ( 0.0 == norm(vecF) )
				msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "Reached vecF = 0.0." );
			break;
		endif
		if ( 0 < iterCount )
			if ( reldiff( vecF, vecF_prev ) == 0.0 )
				msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "There was no change in vecF." );
				break;
			endif
		endif
		%
		iterCount++;
		if ( iterCount > iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "Reached iterMax." );
			break;
		endif
		%
		%
		%
		switch( stepType )
		case STEP_TYPE__BLIND_NEWTON
			vecDelta = - matJ \ vecF;
			vecX_next = vecX + vecDelta;
			vecF_next = funchF( vecX_next ); fevalCount++;
		otherwise
			error( "Invalid value of stepType." );
		endswitch
		%
		%
		%
		switch( jupdateType )
		case JUPDATE_TYPE__NONE
			matJ_next = matJ;
		case JUPDATE_TYPE__BROYDEN
			fooX = vecX_next - vecX;
			fooY = vecF_next - ( vecF + matJ*fooX );
			if ( 0.0 == fooX'*fooX )
				msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "Step size was zero." );
				break;
			endif
			fooJ = (fooY*(fooX'))/(fooX'*fooX);
			matJ_next = matJ + fooJ;
			if ( valdLev >= VALLEV__HIGH )
				assert( reldiff( matJ_next*(vecX_next-vecX), vecF_next-vecF ) <= sqrt(eps) );
			endif
			clear fooX;
			clear fooY;
			clear fooJ;
		case JUPDATE_TYPE__RECALC
			matJ_next = jacobs( vecX, funchF ); jevalCount++;
		otherwise
			error( "Invalid value of jupdateType." );
		endswitch
		%
		%
		%
		vecX_prev = vecX;
		vecF_prev = vecF;
		matJ_prev = matJ;
		%
		vecX = vecX_next;
		vecF = vecF_next;
		matJ = matJ_next;
	endwhile
endfunction

%!function [ vecF, matJ ] = funcFJ_cubyDiagTest( vecX, c )
%!	sizeX = size(vecX,1);
%!	vecXE = (1:sizeX)';
%!	vecF = (vecX-vecXE) + c*(vecX-vecXE).^3;
%!	if ( nargout >= 2 )
%!		matJ = eye(sizeX) + c*3.0*diag((vecX-vecXE).^2);
%!	endif
%!endfunction


%!test
%!	clear;
%!	commondefs;
%!	findLocMin_gnostic_jupdate__defs;
%!	setprngstates(0);
%!	numFigs = 0;
%!	%
%!	caseNum = 30;
%!	msg( __FILE__, __LINE__, sprintf( "caseNum = %d.", caseNum ) );
%!	switch (caseNum)
%!	case 0
%!		sizeX = 2;
%!		c_cuby = 0.0;
%!	case 10
%!		sizeX = 15;
%!		c_cuby = 0.0;
%!	case 20
%!		sizeX = 2;
%!		c_cuby = 0.01;
%!	case 30
%!		sizeX = 15;
%!		c_cuby = 0.01;
%!	case 40
%!		sizeX = 15;
%!		c_cuby = 0.1;
%!	case 100
%!		sizeX = 15;
%!		c_cuby = 0.1;
%!	case 200
%!		sizeX = 20;
%!		c_cuby = 1.0;
%!	otherwise
%!		error( "Ivalid caseNum." );
%!	endswitch
%!	funchFJ = @(dummyX)( funcFJ_cubyDiagTest( dummyX, c_cuby ) );
%!	vecX0 = zeros(sizeX,1);
%!	%
%!	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__NONE ~~~ " );
%!	prm = [];
%!	prm.jupdateType = JUPDATE_TYPE__NONE;
%!	[ vecXF, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchFJ, prm );
%!	vecX0 = zeros(sizeX,1);
%!	%
%!	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__BROYDEN ~~~ " );
%!	prm = [];
%!	prm.jupdateType = JUPDATE_TYPE__BROYDEN;
%!	[ vecXF, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchFJ, prm );
%!	vecX0 = zeros(sizeX,1);
%!	%
%!	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__RECALC ~~~ " );
%!	prm = [];
%!	prm.jupdateType = JUPDATE_TYPE__RECALC;
%!	[ vecXF, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchFJ, prm );
