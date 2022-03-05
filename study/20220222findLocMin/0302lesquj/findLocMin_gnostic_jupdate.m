% Function...
%  A 'gnostic for studying how various forms of Jacobian updating impact convergence.

function [ vecX, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchF, prm=[] )
	commondefs;
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
	%
	STEP_TYPE__BLIND_NEWTON = 0;
	STEP_TYPE__BLIND_GRAD_MIN = 10;
	STEP_TYPE__SCAN_LEV_MIN = 100;
	stepType = mygetfield( prm, "stepType", STEP_TYPE__BLIND_NEWTON );
	if ( valdLev >= VALLEV__LOW )
		assert( isrealscalar(stepType) );
	endif
	%
	JUPDATE_TYPE__NONE = 0;
	JUPDATE_TYPE__BROYDEN = 10;
	JUPDATE_TYPE__REORTHONORM = 20;
	JUPDATE_TYPE__LESQUJ_PRIMAL = 30;
	JUPDATE_TYPE__RECALC = 100;
	jupdateType = mygetfield( prm, "jupdateType", JUPDATE_TYPE__NONE ); % THIS WORKS BETTER???
	jupdateType = mygetfield( prm, "jupdateType", JUPDATE_TYPE__BROYDEN );
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
			msg( __FILE__, __LINE__, sprintf( "  %10.3e, %4d;  %5d,  %3d;  %10.3e.", ...
			  time()-time0, iterCount, ...
			  fevalCount, jevalCount, ...
			  sumsq(vecF)/2.0 ) );
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
			fooJ = (fooY*(fooX'))/(fooX'*fooX);
			matJ_next = matJ + fooJ;
			if ( valdLev >= VALLEV__HIGH )
				assert( reldiff( matJ_next*(vecX_next-vecX), vecF_next-vecF ) <= sqrt(eps) );
			endif
			clear fooX;
			clear fooY;
			clear fooJ;
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
%!	setprngstates(0);
%!	numFigs = 0;
%!	%
%!	sizeX = 15;
%!	c_cuby = 0.1;
%!	funchFJ = @(dummyX)( funcFJ_cubyDiagTest( dummyX, c_cuby ) );
%!	vecX0 = zeros(sizeX,1);
%!	%
%!	prm = [];
%!	[ vecXF, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchFJ, prm );
