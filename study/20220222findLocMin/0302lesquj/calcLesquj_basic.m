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
	useDistanceWeights = true;
	useLatestPtAs0 = true;
	if ( ~isempty(prm) )
		useDistanceWeights = mygetfield( prm, "useDistanceWeights", useDistanceWeights );
	endif
	if ( debugMode )
		assert( isscalar(useDistanceWeights) );
		assert( isbool(useDistanceWeights) );
	endif
	matIX = eye(sizeX,sizeX);
	%
	%
	%
	if (useLatestPtAs0)
		indexOfPt0 = numVals;
	else
		[ foo, indexOfPt0 ] = min( sumsq(vecFVals,1) );
	endif
	vecX0 = vecXVals( :, indexOfPt0 );
	vecF0 = vecFVals( :, indexOfPt0 );
	% Trimming off "0", we'll use the term "pts" instead of "vals".
	vecXPts = [ vecXVals(:,1:indexOfPt0-1), vecXVals(:,indexOfPt0+1:end) ];
	vecFPts = [ vecFVals(:,1:indexOfPt0-1), vecFVals(:,indexOfPt0+1:end) ];
	numPts = size(vecXPts,2);
	%
	%
	%
	vecDeltaPts = vecXPts - vecX0; % sizeX x numPts
	deltaNormSqPts = sumsq(vecDeltaPts,1);
	minDeltaNormSq = min(deltaNormSqPts) + (eps^2)*max(deltaNormSqPts);
	%
	% If we have a jeval, we'll hack it it as a cluster of points.
	% This seems a bit silly, but, it's simple,
	% And I see no obvious reason why it'd be wrong.
	jevalDat = mygetfield( prm, "jevalDat", [] );
	if (~isempty(jevalDat))
		numJevals = size(jevalDat);
		for n=1:numJevals
			vecDeltaX_jeval = minDeltaNormSq*matIX;
			% Octave doesn't auto-broadcast for a "diagonal" matrix.
			% The mathematically pointless "+0.0" gets Octave to convert the matrix to a full matrix.
			% Note that the feval of the jeval is not incorporated here;
			% it should be incorporated as a  separate feval.
			% This is because this way is convenient for me at the moment.
			vecX_jeval = jevalDat(n).vecX + ( vecDeltaX_jeval+0.0 );
			vecF_jeval = jevalDat(n).vecF + ( jevalDat(n).matJ*vecDeltaX_jeval + 0.0 );
			vecXPts = [ vecXPts, vecX_jeval ];
			vecFPts = [ vecFPts, vecF_jeval ];
		endfor
		vecDeltaPts = vecXPts - vecX0; % sizeX x numPts
		deltaNormSqPts = sumsq(vecDeltaPts,1);
		numPts = size(vecXPts,2);
	endif
	%
	%
	%
	if (useDistanceWeights)
		assert( isrealscalar(minDeltaNormSq) );
		assert( 0.0 < minDeltaNormSq );
		assert( isrealarray(deltaNormSqPts) );
		wDistPts = (2.0 * minDeltaNormSq) ./ ( minDeltaNormSq + deltaNormSqPts );
		assert( isrealarray(wDistPts) );
		wResPts = norm(vecF0)./norm(vecFPts);
		wResPts .^= 0.5;
		wDistPts .^= 0.5;
		wPts = wDistPts.*wResPts;
		assert( isrealarray(wPts) );
		vecRhoPts = vecFPts - vecF0;
		foo = vecDeltaPts*diag(wPts);
		matA = foo*(vecDeltaPts');
		matY = foo*(vecRhoPts');
	elseif (0)
		% With ^4, wDistPts may get too small in cases.
		wDistPts = (2.0 * minDeltaNormSq) ./ ( minDeltaNormSq + deltaNormSqPts );
		wDistPts .^= 4;
		vecRhoPts = vecFPts - vecF0;      % sizeX x numPts
		foo = vecDeltaPts*diag(wDistPts); % sizeX x numPts
		matA = foo*(vecDeltaPts');        % sizeX x sizeX
		matY = foo*(vecRhoPts');          % sizeX x sizeF
	else
		vecRhoPts = vecFPts - vecF0;       % sizeX x numPts
		matA = vecDeltaPts*(vecDeltaPts'); % sizeX x sizeX
		matY = vecDeltaPts*(vecRhoPts');   % sizeX x sizeF
	endif
	%
	%
	safeRelTol = eps^0.35;
	[ matR, cholFlag ] = chol( matA );
	if ( 0~=cholFlag || min(diag(matR)) <= safeRelTol*max(abs(diag(matR))) )
		aScale = max(abs(diag(matA)));
		if ( 0 == aScale )
			error( "matA is zero." );
		endif
		[ matR, cholFlag ] = chol( matA + aScale*matIX );
		if ( 0~=cholFlag || min(diag(matR)) <= safeRelTol*max(abs(diag(matR))) )
			error( "Failed to regularize matA." );
		endif
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
%!	funchF = @(dummyX)( vecF_secret + matJ_secret*(dummyX-vecX_secret) + matEta_secret*((dummyX-vecX_secret).^3) );
%!	%
%!	numVals = sizeX+1;
%!	vecXVals = randn(sizeX,numVals);
%!	vecXVals = [ vecXVals, 10*randn(sizeX,2*numVals) ];
%!	vecFVals = funchF(vecXVals) + sigmaLevel*randn(size(vecXVals));
%!	numVals = size(vecXVals,2)
%!	%plot( sumsq(vecFVals,1), 'o-' );
%!	%
%!	prm = [];
%!	prm.jevalDat(1).vecX = randn(sizeX,1);
%!	prm.jevalDat(1).vecF = funchF( prm.jevalDat(1).vecX );
%!	prm.jevalDat(1).matJ = matJ_secret;
%!	[ vecX0, vecF0, matJ0, datOut ] = calcLesquj_basic( vecXVals, vecFVals, prm );
%!	%
%!	vecX0_basic = vecX0;
%!	vecXRoot_basic = vecX0 - matJ0 \ vecF0;
%!	vecXRoot_secret = vecX_secret - matJ_secret \ vecF_secret;
%!	%
%!	vecF0_basic = funchF( vecX0 );
%!	vecFRoot_basic = funchF( vecXRoot_basic );
%!	vecFRoot_secret = funchF( vecXRoot_secret );
%!	%
%!	omega0_basic = sumsq(vecF0_basic)/2.0
%!	omegaRoot_basic = sumsq(vecFRoot_basic)/2.0
%!	omegaRoot_secret = sumsq(vecFRoot_secret)/2.0
