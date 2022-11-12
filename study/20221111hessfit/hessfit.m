function [ f, vecG, matH, datOut ] = hessfit( sizeX, numPts, matX, rvecF, matG, prm=[] )
	%
	[ datCnsts, vecFGL0, hess2lambdaDat ] = __init( sizeX, numPts, matX, rvecF, matG, prm );
	datOut = [];
	%
	%
	[ f, vecG, matH ] = __unpackFGL( sizeX, vecFGL0, hess2lambdaDat );
	msg( __FILE__, __LINE__, "Goodbye." );
	return;
	%
	msg( __FILE__, __LINE__, "TODO: Add options to keep f and even vecG constant." );
	%
	funchRes = @(fgl)( __funcResFGL( fgl ) );
	vecRes0 = funchRes( vecFGL0 );
	funchA = @(fgl)( funchRes(vecFGL0+fgl) - vecRes0 );
	%
	gmres_RTOL = mygetfield( prm, "gmres_RTOL", 1.0e-12 );
	gmres_MAXIT = mygetfield( prm, "gmres_MAXIT", numUnk );
	msg( __FILE__, __LINE__, "Calling gmres()..." );
	%[ vecFGL_delta, gmres_FLAG, gmres_RELRES, gmres_ITER, gmres_RESVEC ] = gmres( funchA, -vecRes0, [], gmres_RTOL, gmres_MAXIT );
	vecFGL_delta = gmres( funchA, -vecRes0, [], gmres_RTOL, gmres_MAXIT )
	vecFGL = vecFGL0 + vecFGL_delta;
	%
	[ f, vecG, matH ] = __unpackFGL( sizeX, vecFGL, hess2lambdaDat );
return;
endfunction


function [ datCnsts, vecFGL0, hess2lambdaDat ] = __init( sizeX, numPts, matX, rvecF, matG, prm )
	assert( 1 <= sizeX );
	assert( 1 <= numPts );
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	%
	datCnsts.ary3XXT = zeros(sizeX,sizeX,numPts);
	for p = 1 : numPts
		datCnsts.ary3XXT(:,:,p) = matX(:,p) * (matX(:,p)');
	endfor
	clear p;
	%
	datCnsts.rvecW0 = mygetfield( prm, "rvecW0", ones(1,numPts) );
	assert( isrealarray(datCnsts.rvecW0,[1,numPts]) );
	%
	datCnsts.rvecW1 = mygetfield( prm, "rvecW1", ones(1,numPts) );
	assert( isrealarray(datCnsts.rvecW1,[1,numPts]) );
	datCnsts.s1 = sum( datCnsts.rvecW1, 2 );
	datCnsts.vecS1X = sum( datCnsts.rvecW1 .* matX, 2 ); % Autobroadcast rvecW1.
	datCnsts.vecS1G = sum( datCnsts.rvecW1 .* matG, 2 ); % Autobroadcast rvecW1.
	datCnsts.matS1XXT = sum( reshape(datCnsts.rvecW1,[1,1,numPts]) .* datCnsts.ary3XXT, 3 ); % Autobroadcast.
	datCnsts.matS1XGTPGXT = zeros( sizeX, sizeX );
	for p = 1 : numPts
		datCnsts.matS1XGTPGXT += ( datCnsts.rvecW1(p) * matX(:,p) ) * ( matG(:,p)' );
	endfor
	clear p;
	datCnsts.matS1XGTPGXT += (datCnsts.matS1XGTPGXT');
	%
	datCnsts.epsHRegu = mygetfield( prm, "epsHRegu", 0.0 );
	assert( isrealscalar(datCnsts.epsHRegu) );
	%
	f0 = mygetfield( prm, "f0", 0.0 );
	vecG0 = mygetfield( prm, "vecG0", zeros(sizeX,1) );
	matH0 = mygetfield( prm, "matH0", zeros(sizeX,sizeX) );
	assert( isrealscalar(f0) );
	assert( isrealarray(vecG0,[sizeX,1]) );
	assert( isrealarray(matH0,[sizeX,sizeX]) );
	assert( issymmetric(matH0) );
	[ vecLambda0, hess2lambdaDat ] = hess2lambda( matH0 );
	vecFGL0 = [ f0; vecG0; vecLambda0 ];
return;
endfunction


function [ f, vecG, matH ] = __unpackFGL( sizeX, vecFGL, hess2lambdaDat )
	f = vecFGL(1);
	vecG = vecFGL(2:1+sizeX);
	matH = lambda2hess( vecFGL(2+sizeX:end), hess2lambdaDat );
return;
endfunction
