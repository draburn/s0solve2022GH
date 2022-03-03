% Function...

function [ vecX0, vecF0, matJ0, datOut ] = calcLesquj_basic( vecXVals, vecFVals, prm=[] )
	%
	sizeX = size(vecXVals,1);
	sizeF = size(vecFVals,1);
	numVals = size(vecXVals,2);
	%
	debugMode = mygetfield( prm, "debugMode", false );
	if ( debugMode )
		msg( __FILE__, __LINE__, "Using debugMode." );
		assert( isrealarray(vecXVals,[sizeX,numVals]) );
		assert( isrealarray(vecFVals,[sizeF,numVals]) );
	endif
	%
	omegaVals = sumsq(vecFVals,1)/2.0;
	[ sorted_omegaVals, sorted_indexVals ] = sort( omegaVals );
	indexOfPt0 = sorted_indexVals(1);
	vecX0 = vecXVals( :, indexOfPt0 );
	vecF0 = vecFVals( :, indexOfPt0 );
	% Trimming off "0", we'll use the term "pts" instead of "vals".
	vecXPts = vecXVals( :, sorted_indexVals(2:end) );
	vecFPts = vecFVals( :, sorted_indexVals(2:end) );
	numPts = numVals - 1;
	clear sorted_omegaVals;
	clear sorted_indexVals;
	%
	vecDeltaPts = vecXPts - vecX0; % sizeX x numPts
	vecRhoPts = vecFPts - vecF0; % sizeX x numPts
	matA = vecDeltaPts*(vecDeltaPts'); % sizeX x sizeX
	matY = vecDeltaPts*(vecRhoPts'); % sizeX x sizeF
	%
	[ matR, cholFlag ] = chol( matA );
	if ( 0~=cholFlag )
		aScale = max(abs(diag(matA)));
		matR = chol( matA + aScale*eye(sizeX,sizeX) );
	endif
	matJ0T = matR \ ( matR' \ matY );
	matJ0 = matJ0T';
	if ( nargout >= 2 )
		datOut.matRes = ( matY - (matA*matJ0T) )';
	endif
return;
endfunction


%!test
%!	setprngstates(0);
%!	sizeX = 5;
%!	sizeF = 5;
%!	vecX_secret = randn(sizeX,1)
%!	vecF_secret = randn(sizeF,1)
%!	matJ_secret = randn(sizeF,sizeX)
%!	funchF = @(dummyX)( vecF_secret + matJ_secret*(dummyX-vecX_secret) );
%!	numVals = 100;
%!	vecXVals = randn(sizeX,numVals);
%!	vecFVals = funchF(vecXVals);
%!	%
%!	prm = [];
%!	[ vecX0, vecF0, matJ0, datOut ] = calcLesquj_basic( vecXVals, vecFVals, prm )
