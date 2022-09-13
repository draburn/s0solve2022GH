% Function...

function [ funchFOfX, datOut ] = genFunchAPS2022( bigN, prm=[] )
	time_start = time();
	% bigN: num unk, x-size, domain size, num col in J, etc.
	assert( isposintscalar(bigN) );
	%
	% bigM: f-size, range size, "measure" size, num row in J, etc.
	bigM = bigN;
	bigM = mygetfield( prm, "bigM", bigM );
	assert( isposintscalar(bigM) );
	%
	% lambda: average num "Large Elements" per row.
	lambda = sqrt(bigN);
	lambda = mygetfield( prm, "lambda", lambda );
	assert( isrealscalar(lambda) );
	assert( 0.0 == lambda || 1.0 <= lambda );
	%
	% bigL: total number of "Large Elements".
	% Require either none or at least one per row.
	bigL = ceil( lambda * bigN );
	assert( 0 == bigL || bigN <= bigL );
	%
	% bigP: "compressed" size for MxN matrices.
	bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
	bigP = mygetfield( prm, "bigP", bigP );
	assert( isposintscalar(bigP) );
	assert( bigP <= bigM );
	assert( bigP <= bigN );
	%
	% Check anticipated memory usage.
	% Note that inclusion of data in datOut is free.
	% But, about 50% additional memory is required during construction.
	allowedMemInGB = 12.0;
	allowedMemInGB = mygetfield( prm, "allowedMemInGB", allowedMemInGB );
	assert( isrealscalar(allowedMemInGB) );
	assert( 0.0 < allowedMemInGB );
	numMxPxPxN = 7;
	bytesPerDouble = 8;
	bytesPerMxPxPxN = bytesPerDouble * ( bigM*bigP + bigP*bigN );
	numSparseMxN = 6;
	bytesPerSparseElem = 16.0; % Approximate.
	bytesPerSparseMxN = ceil( bytesPerSparseElem * bigL );
	bytesRequired = numMxPxPxN*bytesPerMxPxPxN + numSparseMxN*bytesPerSparseMxN;
	estMemInGB = bytesRequired / 1.0E9;
	if ( estMemInGB > allowedMemInGB )
		error([ "Required memory (" num2str(estMemInGB) "GB) exceeds allowed memory (" num2str(allowedMemInGB) "GB)." ]);
	endif
	%
	%
	% Get coefficients.
	cx = mygetfield( prm, "cx", 1.0E0 );
	c0 = mygetfield( prm, "c0", 0.0E0 );
	c1 = mygetfield( prm, "c1", 1.0E0 );
	c2 = mygetfield( prm, "c2", 0.0E0 );
	c3 = mygetfield( prm, "c3", 0.0E0 );
	s1 = mygetfield( prm, "s1", 0.0E0 );
	s2 = mygetfield( prm, "s2", 0.0E0 );
	s3 = mygetfield( prm, "s3", 0.0E0 );
	assert( isrealscalar(cx) );
	assert( isrealscalar(c0) );
	assert( isrealscalar(c1) );
	assert( isrealscalar(c2) );
	assert( isrealscalar(c3) );
	assert( isrealscalar(s1) );
	assert( isrealscalar(s2) );
	assert( isrealscalar(s3) );
	assert( 0.0 < cx );
	assert( 0.0 <= c0 );
	assert( 0.0 <= c1 );
	assert( 0.0 <= c2 );
	assert( 0.0 <= c3 );
	assert( 0.0 <= s1 );
	assert( 0.0 <= s2 );
	assert( 0.0 <= s3 );
	%
	%
	% Generate values.
	prngStateDat = getprngstatedat();
	rSeed = mygetfield( prm, "rSeed", [] );
	if ( ~isempty(rSeed) )
		setprngstates(rSeed);
	end
	time_startBuild = time();
	%
	vecXSecret = randn( bigN, 1 )*cx;
	%
	matA11L = randn( bigM, bigP )*c1;
	matA11R = randn( bigP, bigN );
	matA12L = randn( bigM, bigP )*c1;
	matA12R = randn( bigP, bigN );
	matA21L = randn( bigM, bigP )*c2;
	matA21R = randn( bigP, bigN );
	matA22L = randn( bigM, bigP )*c2;
	matA22R = randn( bigP, bigN );
	matA31L = randn( bigM, bigP )*c3;
	matA31R = randn( bigP, bigN );
	matA32L = randn( bigM, bigP )*c3;
	matA32R = randn( bigP, bigN );
	matA33L = randn( bigM, bigP )*c3;
	matA33R = randn( bigP, bigN );
	%
	if ( 0 == bigL )
		vecM = [];
		vecN = [];
		matB1  = sparse( [], [], [], bigM, bigN );
		matB21 = sparse( [], [], [], bigM, bigN );
		matB22 = sparse( [], [], [], bigM, bigN );
		matB31 = sparse( [], [], [], bigM, bigN );
		matB32 = sparse( [], [], [], bigM, bigN );
		matB33 = sparse( [], [], [], bigM, bigN );
	else
		assert( bigL >= bigM ); % Require at least one per row.
		vecM = 1 + floor( (bigM-eps)*rand(bigL,1) );
		vecM(1:bigM) = (1:bigM)'; % Force at least one per row.
		vecN = 1 + floor( (bigN-eps)*rand(bigL,1) );
		% Note that some elements may be repeated, so there may be fewer than bigL "Large Entries".
		matB1  = sparse( vecM, vecN, s1*randn(bigL,1), bigM, bigN );
		matB21 = sparse( vecM, vecN, s2*randn(bigL,1), bigM, bigN );
		matB22 = sparse( vecM, vecN, s2*randn(bigL,1), bigM, bigN );
		matB31 = sparse( vecM, vecN, s3*randn(bigL,1), bigM, bigN );
		matB32 = sparse( vecM, vecN, s3*randn(bigL,1), bigM, bigN );
		matB33 = sparse( vecM, vecN, s3*randn(bigL,1), bigM, bigN );
	endif
	%
	sqrtM = sqrt(bigM);
	funchYOfX = @(x)(  (x - vecXSecret) / sqrtM );
	funchF1OfY = @(y)( ...
	   ( matA11L*(matA11R*y) ) + ( matA12L*(matA12R*y) ) ...
	 + ( (matA21L*(matA21R*y)) .* (matA22L*(matA22R*y)) ) ...
	 + ( (matA31L*(matA31R*y)) .* (matA32L*(matA32R*y)) .* (matA33L*(matA33R*y)) ) ...
	 + ( matB1*y ) ...
	 + ( (matB21*y) .* (matB22*y) ) ...
	 + ( (matB31*y) .* (matB32*y) .* (matB33*y) ) ...
	);
	
	if ( bigM == bigN )
		funchFOfX = @(x)(  funchF1OfY(funchYOfX(x)) + (c0*funchYOfX(x))  );
	else
		vecInterpN = (0:bigN-1)/(bigN-1.0);
		vecInterpM = (0:bigM-1)/(bigM-1.0);
		funchYInterpOfX = @(x)(  interp1( vecInterpN, funchYOfX(x), vecInterpM', "linear" )  );
		funchFOfX = @(x)(  funchF1OfY(funchYOfX(x)) + (c0*funchYInterpOfX(x))  );
		% It'd be nice to have a matrix of this, but, meh.
	endif
	%
	time_finishBuild = time();
	setprngstatedat(prngStateDat);
	%
	%
	%
	if (nargout>=2)
		% The basics.
		datOut.bigN = bigN;
		datOut.bigM = bigM;
		% The big secret.
		datOut.vecXSecret = vecXSecret;
		% Locations of "Large Elements".
		datOut.vecMSparse = vecM;
		datOut.vecNSparse = vecN;
		% Sufficient data to construct Jacobian near root.
		datOut.c0 = c0;
		datOut.matA11L = matA11L;
		datOut.matA11R = matA11R;
		datOut.matA12L = matA12L;
		datOut.matA12R = matA12R;
		datOut.matB1 = matB1;
		% Misc.
		datOut.estMemInGB = estMemInGB;
		datOut.initTimeInS = time_startBuild - time_start;
		datOut.buildTimeInS = time_finishBuild - time_startBuild;
	endif
return;
endfunction
