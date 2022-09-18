function [ vecXBest, grootFlag, fevalCount, datOut ] = groot_findZero_800( funchF, vecX0, prm=[] )
	if ( 0 == nargin )
		vecXBest = __FILE__;
		return;
	elseif ( nargin < 2 )
		error( "Too few input arguments." );
	elseif ( 3 < nargin )
		error( "Too many input arguments." );
	elseif ( 4 < nargout )
		error( "Too many output arguments." );
	endif
	groot__commonInit;
	vecXBest = [];
	grootFlag = GROOT_FLAG__VOID;
	fevalCount = 0;
	datOut = [];
	%
	%
	findZero800Prm = mygetfield( prm, "findZero800Prm", [] );
	findZero800Prm.verbLev = mygetfield( findZero800Prm, "verbLev", VERBLEV__FLAGGED );
	[ vecXBest, vecFF, findZero800DatOut ] = findZero_800( vecX0, funchF, findZero800Prm );
	fevalCount = findZero800DatOut.fevalCount;
	if (  fevalCount <= prm.fevalLimit  &&  norm(vecFF) <= prm.fTol*(1.0+100.0*eps)  )
		grootFlag = GROOT_FLAG__CNVG;
	elseif ( fevalCount >= prm.fevalLimit )
		grootFlag = GROOT_FLAG__STOP;
	else
		grootFlag = GROOT_FLAG__FAIL;
	endif
	%
	datOut.elapsedTime = 0.0;
	datOut.matRecordX = [];
	datOut.matInfoA = [];
	datOut.matInfoB = [];
return;
endfunction