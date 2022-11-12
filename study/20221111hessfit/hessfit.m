function [ f, vecG, matH, datOut ] = hessfit( sizeX, numPts, matX, rvecF, matG, prm=[] )
	%
	[ datCnsts, vecFGL0, hess2lambdaDat ] = __init( sizeX, numPts, matX, rvecF, matG, prm );
	datOut = [];
	%
	rTol = mygetfield( prm, "rTol", 1.0e-12 );
	useCnstF = mygetfield( prm, "useCnstF", false );
	useCnstG = mygetfield( prm, "useCnstG", false ); % subsumes useCnstF.
	if ( useCnstG && useCnstF )
		vecFG0 = vecFGL0(1:1+sizeX);
		vecL0 = vecFGL0(2+sizeX:end);
		funchRes = @(vecL)( __funcResFGL( [vecFG0;vecL], sizeX, numPts, matX, rvecF, matG, datCnsts, hess2lambdaDat )(2+sizeX:end) );
		vecRes0 = funchRes( vecL0 );
		funchA = @(vecL)( funchRes(vecL0+vecL) - vecRes0 );
		vecB = -vecRes0;
		numUnk = (sizeX*(sizeX+1))/2;
		%
		maxIt = mygetfield( prm, "maxIt", numUnk  );
		assert( 0.0 < rTol );
		assert( isposintscalar(maxIt) );
		vecGL_delta = cgs( funchA, vecB, rTol, maxIt );
		%
		vecFGL = vecFGL0;
		vecFGL(2+sizeX:end) += vecGL_delta;
	elseif ( useCnstF )
		f0 = vecFGL0(1);
		vecGL0 = vecFGL0(2:end);
		funchRes = @(gl)( __funcResFGL( [f0;gl], sizeX, numPts, matX, rvecF, matG, datCnsts, hess2lambdaDat )(2:end) );
		vecRes0 = funchRes( vecGL0 );
		funchA = @(gl)( funchRes(vecGL0+gl) - vecRes0 );
		vecB = -vecRes0;
		numUnk = (sizeX*(sizeX+1))/2 + sizeX;
		%
		maxIt = mygetfield( prm, "maxIt", numUnk  );
		assert( 0.0 < rTol );
		assert( isposintscalar(maxIt) );
		vecGL_delta = cgs( funchA, vecB, rTol, maxIt );
		%
		vecFGL = vecFGL0;
		vecFGL(2:end) += vecGL_delta;
	elseif ( useCnstG )
		error( "Use of useCnstG without useCnstF is not supported." );
	else
		funchRes = @(fgl)( __funcResFGL( fgl, sizeX, numPts, matX, rvecF, matG, datCnsts, hess2lambdaDat ) );
		vecRes0 = funchRes( vecFGL0 );
		funchA = @(fgl)( funchRes(vecFGL0+fgl) - vecRes0 );
		vecB = -vecRes0;
		numUnk = (sizeX*(sizeX+1))/2 + sizeX + 1;
		%
		maxIt = mygetfield( prm, "maxIt", numUnk  );
		assert( 0.0 < rTol );
		assert( isposintscalar(maxIt) );
		vecFGL_delta = cgs( funchA, vecB, rTol, maxIt );
		%
		vecFGL = vecFGL0 + vecFGL_delta;
	endif
	[ f, vecG, matH ] = __unpackFGL( sizeX, vecFGL, hess2lambdaDat );
	%
	if ( nargout >= 4 )
		[ f, vecG, matH ] = __unpackFGL( sizeX, vecFGL, hess2lambdaDat );
		%
		datOut.rvecRho =  f + (vecG'*matX) + (sum( matX .* ( matH * matX ), 1 )/2.0) - rvecF; % Autobroadcast vecG.
		datOut.rvecW0Rho = datCnsts.rvecW0 .* datOut.rvecRho;
		datOut.matR = vecG + ( matH * matX ) - matG;
		datOut.matW1R = datCnsts.rvecW1 .* datOut.matR; % Autobroadcast rvecW1.
		datOut.rvecRNorm = sqrt( sum(datOut.matR.^2,1) );
		datOut.rvecW1RNorm = sqrt( sum(datOut.matW1R.^2,1) );
		%
		datOut.stats.rho = [ mean(datOut.rvecRho), std(datOut.rvecRho) ];
		datOut.stats.absRho = [ mean(abs(datOut.rvecRho)), std(abs(datOut.rvecRho)) ];
		datOut.stats.rNorm = [ mean(datOut.rvecRNorm), std(datOut.rvecRNorm) ];
		datOut.stats.hFrobNorm = sqrt(sum(sum(matH.^2)));
		datOut.stats.w0rho = [ mean(datOut.rvecW0Rho), std(datOut.rvecW0Rho) ];
		datOut.stats.absW0Rho = [ mean(abs(datOut.rvecW0Rho)), std(abs(datOut.rvecW0Rho)) ];
		datOut.stats.w1RNorm = [ mean(datOut.rvecW1RNorm), std(datOut.rvecW1RNorm) ];
		datOut.stats.hReguFrobNorm = datCnsts.epsHRegu * datOut.stats.hFrobNorm;
	endif
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
	datCnsts.rvecW1 = mygetfield( prm, "rvecW1", ones(1,numPts)*(eps^0.4)/sizeX );
	%%%datCnsts.rvecW1 = mygetfield( prm, "rvecW1", ones(1,numPts) );
	assert( isrealarray(datCnsts.rvecW1,[1,numPts]) );
	datCnsts.s1 = sum( datCnsts.rvecW1, 2 );
	datCnsts.vecS1X = sum( datCnsts.rvecW1 .* matX, 2 ); % Autobroadcast rvecW1.
	datCnsts.vecS1G = sum( datCnsts.rvecW1 .* matG, 2 ); % Autobroadcast rvecW1.
	datCnsts.matS1XXT = sum( reshape(datCnsts.rvecW1,[1,1,numPts]) .* datCnsts.ary3XXT, 3 ); % Autobroadcast thingy.
	datCnsts.matS1XGTPGXT = zeros( sizeX, sizeX );
	for p = 1 : numPts
		datCnsts.matS1XGTPGXT += ( datCnsts.rvecW1(p) * matX(:,p) ) * ( matG(:,p)' );
	endfor
	clear p;
	datCnsts.matS1XGTPGXT += (datCnsts.matS1XGTPGXT');
	%
	datCnsts.epsHRegu = mygetfield( prm, "epsHRegu", (eps^0.8)/(sizeX^2) );
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


function vecResFGL = __funcResFGL( vecFGL, sizeX, numPts, matX, rvecF, matG, datCnsts, hess2lambdaDat )
	[ f, vecG, matH ] = __unpackFGL( sizeX, vecFGL, hess2lambdaDat );
	%
	rvecRho =  f + (vecG'*matX) + (sum( matX .* ( matH * matX ), 1 )/2.0) - rvecF; % Autobroadcast vecG.
	rvecW0Rho = datCnsts.rvecW0 .* rvecRho;
	%
	partialF = sum( rvecW0Rho, 2 );
	vecNablaG = sum( rvecW0Rho .* matX, 2 ) + ( datCnsts.s1 * vecG ) + ( matH * datCnsts.vecS1X ) - datCnsts.vecS1G;
	foo = ( vecG * datCnsts.vecS1X' ) + ( matH * datCnsts.matS1XXT );
	matNablaL = sum( reshape(rvecW0Rho,[1,1,numPts]) .* datCnsts.ary3XXT, 3 ) ... % Autobroadcast thingy.
	  + foo + (foo') + ( datCnsts.epsHRegu * matH ) - datCnsts.matS1XGTPGXT;
	%
	% Oooooo...
	%  should this be vech(matNablaL) or hess2lambda(matNablaL)?
	% It shouldn't really matter, because the root is the same either way.
	vecResFGL = [ partialF; vecNablaG; vech(matNablaL) ];
return;
endfunction
