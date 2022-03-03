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
	if (0)
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
	else
	[ foo, indexOfPt0 ] = min( sumsq(vecFVals,1) );
	vecX0 = vecXVals( :, indexOfPt0 );
	vecF0 = vecFVals( :, indexOfPt0 );
	% Trimming off "0", we'll use the term "pts" instead of "vals".
	vecXPts = [ vecXVals(:,1:indexOfPt0-1), vecXVals(:,indexOfPt0+1:end) ];
	vecFPts = [ vecFVals(:,1:indexOfPt0-1), vecFVals(:,indexOfPt0+1:end) ];
	numPts = numVals - 1;
	endif
	%
	vecDeltaPts = vecXPts - vecX0; % sizeX x numPts
	deltaNormSqPts = sumsq(vecDeltaPts,1);
	minDeltaNormSq = min(deltaNormSqPts);
	wDistPts = (2.0 * minDeltaNormSq) ./ ( minDeltaNormSq + deltaNormSqPts );
	%
	vecRhoPts = vecFPts - vecF0; % sizeX x numPts
	if (0)
	matA = vecDeltaPts*(vecDeltaPts'); % sizeX x sizeX
	matY = vecDeltaPts*(vecRhoPts'); % sizeX x sizeF
	else
	foo = vecDeltaPts*diag(wDistPts); % sizeX x numPts
	matA = foo*(vecDeltaPts');        % sizeX x sizeX
	matY = foo*(vecRhoPts');          % sizeX x sizeF
	endif
	%
	[ matR, cholFlag ] = chol( matA );
	if ( 0~=cholFlag )
		aScale = max(abs(diag(matA)));
		matR = chol( matA + aScale*eye(sizeX,sizeX) );
	endif
	matJ0T = matR \ ( matR' \ matY );
	matJ0 = matJ0T';
	if ( nargout >= 2 )
		%vecResPts = vecRhoPts - matJ0*vecDeltaPts;
		datOut = [];
	endif
return;
endfunction


%!test
%!	setprngstates(0);
%!	sizeX = 20
%!	sizeF = 20
%!	etaLevel = 1e-2
%!	sigmaLevel = 1e-2
%!	%
%!	vecX_secret = randn(sizeX,1);
%!	vecF_secret = 0*randn(sizeF,1);
%!	matJ_secret = randn(sizeF,sizeX);
%!	matEta_secret = etaLevel*randn(sizeF,sizeX);
%!	funchF = @(dummyX)( vecF_secret + matJ_secret*(dummyX-vecX_secret) + matEta_secret*((dummyX-vecX_secret).^2) );
%!	%
%!	numVals = sizeX+2
%!	vecXVals = randn(sizeX,numVals);
%!	vecFVals = funchF(vecXVals) + sigmaLevel*randn(size(vecXVals));
%!	%
%!	prm = [];
%!	[ vecX0, vecF0, matJ0, datOut ] = calcLesquj_basic( vecXVals, vecFVals, prm );
%!	%
%!	vecX0_basic = vecX0;
%!	vecXRoot_basic = vecX0 - matJ0 \ vecF0;
%!	vecXRoot_secret = vecX_secret - matJ_secret \ vecF_secret;
%!	%
%!	vecF0_basic = funchF(vecF0);
%!	vecFRoot_basic = funchF( vecXRoot_basic );
%!	vecFRoot_secret = funchF( vecXRoot_secret );
%!	%
%!	omega0_basic = sumsq(vecF0_basic)/2.0
%!	omegaRoot_basic = sumsq(vecFRoot_basic)/2.0
%!	omegaRoot_secret = sumsq(vecFRoot_secret)/2.0
