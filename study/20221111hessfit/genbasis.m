% matV = genbasis( matX, prm=[] )
%
% Generate a reasonable basis subspace about matX(:,1).
%
% DRaburn 2022-11-11.

function matV = genbasis( matX, prm=[] )
	sizeX = size(matX,1);
	numPts = size(matX,2);
	%assert( isposintscalar(sizeX) );
	%assert( isposintscalar(numPts) );
	assert( 2 <= numPts );
	assert( isrealarray(matX,[sizeX,numPts]) );
	%assert( isrealarray(rvecF,[1,numPts]) );
	%assert( isrealarray(matG,[sizeX,numPts]) );
	%
	vecX0 = matX(:,1);
	matDX = matX - vecX0; % Autobroadcast vecX0.
	%
	matU = zeros( sizeX, numPts - 1 );
	for p = 1 : numPts-1
		vecDX = matDX(:,p+1);
		s = norm(vecDX);
		assert( s > 0.0 );
		matU(:,p) = vecDX / s;
	endfor
	%
	looseyness = 0;
	switch (looseyness)
	case 3
		bigKMax = mygetfield( prm, "bigKMax", sizeX );
		uNormSqThresh = mygetfield( prm, "uNormSqThresh", 0.0 );
		uFallSqThresh = mygetfield( prm, "uFallSqThreshs", 0.0 );
	case 2
		bigKMax = mygetfield( prm, "bigKMax", min( 2*(numPts-1), sizeX ) );
		uNormSqThresh = mygetfield( prm, "uNormSqThresh", 0.001^2 );
		uFallSqThresh = mygetfield( prm, "uFallSqThreshs", 0.01^2 );
	case 1
		bigKMax = mygetfield( prm, "bigKMax", min( floor((numPts-1)/2.0), sizeX ) );
		uNormSqThresh = mygetfield( prm, "uNormSqThresh", 0.01^2 );
		uFallSqThresh = mygetfield( prm, "uFallSqThreshs", 0.4^2 );
	case 0
		bigKMax = mygetfield( prm, "bigKMax", min( floor(sqrt((numPts-1)/2.0)), sizeX ) );
		uNormSqThresh = mygetfield( prm, "uNormSqThresh", 0.1^2 );
		uFallSqThresh = mygetfield( prm, "uFallSqThreshs", 0.5^2 );
	otherwise
		error( "Invalid case." );
	endswitch
	assert( isposintscalar(bigKMax) );
	assert( bigKMax <= sizeX );
	assert( isrealscalar(uNormSqThresh) );
	assert( isrealscalar(uFallSqThresh) );
	assert( 0.0 <= uNormSqThresh );
	assert( 0.0 <= uFallSqThresh );
	%
	bigK = 0; % K = size(matV,2).
	matV = zeros( sizeX, 0 );
	for p = 1 : numPts-1
		if ( bigK >= bigKMax )
			break;
		endif
		%
		vecU = matU(:,p);
		uNormSq = sumsq(vecU);
		if ( uNormSq < uNormSqThresh )
			continue;
		endif
		%
		vecVWouldBe = vecU / sqrt(uNormSq);
		matUWouldBe = matU;
		uFallSq = uNormSq;
		for q = p+1 : numPts-1
			s = vecVWouldBe' * matUWouldBe(:,q);
			uFallSq += s^2;
			matUWouldBe(:,q) -= s * vecVWouldBe;
		endfor
		if ( uFallSq < uFallSqThresh )
			continue;
		endif
		%
		bigK++;
		matV = [ matV, vecVWouldBe ];
		matU = matUWouldBe;
		clear vecVWouldBe;
		clear matUWouldBe;
	endfor
	%
	%(matV'*matV+1E8)-1E8
	assert( reldiff(matV'*matV,eye(bigK)) < sqrt(eps) );
return;
endfunction
