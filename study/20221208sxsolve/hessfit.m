% function [ f, vecG, matH, datOut ] = hessfit( matX, rvecF, matG, prm=[] )
%
% Generate a Hessian/quadratic fit to the given values about the origin.
% Intended for cases where there is large noise in f and g.
%
% DRaburn 2022-12-02.


function [ f, vecG, matH, datOut ] = hessfit( matX, rvecF, matG, prm=[] )
	datOut = [];
	%
	[ vecFGL0, precalcDat, hess2lambdaDat ] = __init( matX, rvecF, matG, prm );
	%
	funchRes = @(fgl)( __funcResFGL( fgl, matX, rvecF, matG, precalcDat, hess2lambdaDat ) );
	[ vecFGL, datOut.solveDat ] = __solve( funchRes, vecFGL0, prm );
	%
	sizeX = size(matX,1);
	[ f, vecG, matH ] = __unpackFGL( sizeX, vecFGL, hess2lambdaDat );
return
endfunction



function [ vecFGL0, precalcDat, hess2lambdaDat ] = __init( matX, rvecF, matG, prm )
	sizeX = size(matX,1);
	numPts = size(matX,2);
	assert( isposintscalar(sizeX) );
	assert( isposintscalar(numPts) );
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	%
	if ( numPts <= sizeX )
		msg( __FILE__, __LINE__, "Warning: numPts <= sizeX; result may be ill-behaved." );
	endif
	%
	%
	precalcDat.ary3XXT = zeros(sizeX,sizeX,numPts);
	for p = 1 : numPts
		precalcDat.ary3XXT(:,:,p) = matX(:,p) * (matX(:,p)');
	endfor
	clear p;
	%
	%
	rvecG = sqrt(sum(matG.^2,1));
	rvecW0 = 1.0./( abs(rvecF) + sqrt(eps)*max(abs(rvecF)) );
	rvecW1 = 1.0./( abs(rvecG) + sqrt(eps)*max(abs(rvecG)) );
	epsHRegu = (eps^0.8)*sum(rvecF.^2)/sum(sum(matG.^2));
	precalcDat.rvecW0 = mygetfield( prm, "rvecW0", rvecW0 );
	precalcDat.rvecW1 = mygetfield( prm, "rvecW1", rvecW1 );
	precalcDat.epsHRegu = mygetfield( prm, "epsHRegu", epsHRegu );
	assert( isrealarray(precalcDat.rvecW0,[1,numPts]) );
	assert( isrealarray(precalcDat.rvecW1,[1,numPts]) );
	assert( isrealscalar(precalcDat.epsHRegu) );
	precalcDat.s1 = sum( precalcDat.rvecW1, 2 );
	precalcDat.vecS1X = sum( precalcDat.rvecW1 .* matX, 2 ); % Autobroadcast rvecW1.
	precalcDat.vecS1G = sum( precalcDat.rvecW1 .* matG, 2 ); % Autobroadcast rvecW1.
	precalcDat.matS1XXT = sum( reshape(precalcDat.rvecW1,[1,1,numPts]) .* precalcDat.ary3XXT, 3 ); % Autobroadcast thingy.
	precalcDat.matS1XGTPGXT = zeros( sizeX, sizeX );
	for p = 1 : numPts
		precalcDat.matS1XGTPGXT += ( precalcDat.rvecW1(p) * matX(:,p) ) * ( matG(:,p)' );
	endfor
	clear p;
	precalcDat.matS1XGTPGXT += (precalcDat.matS1XGTPGXT');
	%
	%
	f0 = mygetfield( prm, "f0", 0.0 );
	vecG0 = mygetfield( prm, "vecG0", zeros(sizeX,1) );
	matH0 = mygetfield( prm, "matH0", zeros(sizeX,sizeX) );
	assert( isrealscalar(f0) );
	assert( isrealarray(vecG0,[sizeX,1]) );
	assert( isrealarray(matH0,[sizeX,sizeX]) );
	assert( issymmetric(matH0) );
	[ vecLambda0, hess2lambdaDat ] = __hess2lambda( matH0 );
	vecFGL0 = [ f0; vecG0; vecLambda0 ];
return
endfunction



function vecResFGL = __funcResFGL( vecFGL, matX, rvecF, matG, precalcDat, hess2lambdaDat )
	%msg( __FILE__, __LINE__, "Greetings from __funcResFGL()!" );
	sizeX = size(matX,1);
	numPts = size(matX,2);
	[ f, vecG, matH ] = __unpackFGL( sizeX, vecFGL, hess2lambdaDat );
	%
	rvecRho =  f + (vecG'*matX) + (sum( matX .* ( matH * matX ), 1 )/2.0) - rvecF; % Autobroadcast vecG.
	rvecW0Rho = precalcDat.rvecW0 .* rvecRho;
	%
	partialF = sum( rvecW0Rho, 2 );
	vecNablaG = sum( rvecW0Rho .* matX, 2 ) + ( precalcDat.s1 * vecG ) + ( matH * precalcDat.vecS1X ) - precalcDat.vecS1G;
	foo = ( vecG * precalcDat.vecS1X' ) + ( matH * precalcDat.matS1XXT );
	matNablaL = sum( reshape(rvecW0Rho,[1,1,numPts]) .* precalcDat.ary3XXT, 3 ) ... % Autobroadcast thingy.
	  + foo + (foo') + ( precalcDat.epsHRegu * matH ) - precalcDat.matS1XGTPGXT;
	%
	vecResFGL = [ partialF; vecNablaG; vech(matNablaL) ];
return;
endfunction



function [ vecFGL, solveDat ] = __solve( funchRes, vecFGL0, prm )
	startTime = time();
	%[ vecFGL, solveDat ] = __solve_gmres( funchRes, vecFGL0, prm );
	[ vecFGL, solveDat ] = __solve_fullMat( funchRes, vecFGL0, prm );
	%[ vecFGL, solveDat ] = __solve_sparseMat( funchRes, vecFGL0, prm );
	return;
	vecRes0 = funchRes(vecFGL0);
	funchA = @(fgl)( funchRes(fgl+vecFGL0) - vecRes0 );
	linsolfPrm = [];
	linsolfPrm.tol = mygetfield( prm, "tol", 1e-14 );
	vecDeltaFGL = linsolf( funchA, -vecRes0, zeros(length(vecFGL0),1), linsolfPrm );
	vecFGL = vecDeltaFGL + vecFGL0;
	solveDat.elapsedTime = time()-startTime;
	msg( __FILE__, __LINE__, sprintf( "__solve() time: %f", time()-startTime ) );
return
endfunction
function [ vecFGL, solveDat ] = __solve_gmres( funchRes, vecFGL0, prm )
	vecRes0 = funchRes(vecFGL0);
	funchA = @(fgl)( funchRes(fgl+vecFGL0) - vecRes0 );
	rTol = mygetfield( prm, "rTol", 1e-16 );
	maxIt = mygetfield( prm, "maxIt", length(vecRes0) );
	assert( isrealscalar(rTol) );
	assert( 0.0 < rTol );
	assert( rTol < 1.0 );
	assert( isposintscalar(maxIt) );
	assert( maxIt <= length(vecRes0) );
	[ vecDeltaFGL, statFlag, relres, iterCount, vecRes ] = gmres( funchA, -vecRes0, [], rTol, maxIt );
	vecFGL = vecDeltaFGL + vecFGL0;
	solveDat.statFlag = statFlag;
	solveDat.relres = relres;
	solveDat.iterCount = iterCount;
	solveDat.vecRes = vecRes;
return
endfunction
function [ vecFGL, solveDat ] = __solve_fullMat( funchRes, vecFGL0, prm )
	vecRes0 = funchRes(vecFGL0);
	funchA = @(fgl)( funchRes(fgl+vecFGL0) - vecRes0 );
	sz = length(vecRes0);
	% May take way to much memory.
	matM = zeros(sz,sz);
	for n = 1:sz
		fgl = zeros(sz,1);
		fgl(n) = 1.0;
		matM(:,n) = funchA(fgl);
	endfor
	clear fgl;
	%msg( __FILE__, __LINE__, "About our matrix..." );
	%msg( __FILE__, __LINE__, sprintf("  sz = %d", sz ) );
	%msg( __FILE__, __LINE__, sprintf("  cond() = %0.3e", cond(matM) ) )
	vecFGL = vecFGL0 - matM \ vecRes0;
	solveDat = [];
return
endfunction
function [ vecFGL, solveDat ] = __solve_sparseMat( funchRes, vecFGL0, prm )
	vecRes0 = funchRes(vecFGL0);
	funchA = @(fgl)( funchRes(fgl+vecFGL0) - vecRes0 );
	sz = length(vecRes0);
	matM = spalloc( sz, sz, ceil(100 + 10*sz) );
	% The size is just a guess. nnz seems to be about 4.8 ~ 5.0 * sz.
	for n = 1:sz
		fgl = zeros(sz,1);
		fgl(n) = 1.0;
		vecT = funchA(fgl);
		t = sqrt(sum(vecT.^2)/sz);
		matM(:,n) = vecT.*(abs(vecT)>(eps^0.8)*t);
		% This thresholding might eliminate the regularization for H. Oh well.
	endfor
	clear fgl;
	%msg( __FILE__, __LINE__, "About our matrix..." );
	%msg( __FILE__, __LINE__, sprintf("  sz = %d", sz ) );
	%msg( __FILE__, __LINE__, sprintf("  nnz = %d (%f%%)", nnz(matM), nnz(matM)*100.0/(sz*sz) ) );
	%msg( __FILE__, __LINE__, sprintf("  cond() = %0.3e", cond(matM) ) );
	%full(matM)
	vecFGL = vecFGL0 - matM \ vecRes0;
	solveDat = [];
return
endfunction



function [ vecLambda, datOut ] = __hess2lambda( matHess )
	vecLambda = vech( matHess - diag(diag(matHess)/2.0) );
	if ( 1 == nargout )
		return;
	endif
	sz = size(matHess,1);
	datOut.sz = sz;
	datOut.matDuplish = sparse(duplication_matrix(sz));
	% We need to halve the coeff corresponding to the diagonal of the hessian.
	for n = 1 : sz
		m = (n-1)*sz + n;
		datOut.matDuplish( m, m - (n*(n-1))/2 ) = 2.0;
	endfor
return;
endfunction

function matH = __lambda2hess( vecLambda, datIn )
	matH = reshape( datIn.matDuplish * vecLambda, [datIn.sz,datIn.sz] );
return;
endfunction

function [ f, vecG, matH ] = __unpackFGL( sizeX, vecFGL, hess2lambdaDat )
	% Note: __unpackFGL returns f, g, and *H* from {f, g, lambda}.
	f = vecFGL(1);
	vecG = vecFGL(2:1+sizeX);
	matH = __lambda2hess( vecFGL(2+sizeX:end), hess2lambdaDat );
return;
endfunction


%!test
%!	tic();
%!	%setprngstates(55436400);
%!	setprngstates();
%!	gNoiseLevel = 0.0
%!	fNoiseLevel = 0.0
%!	expVarCoeff = 1.0
%!	forceOrthogonalBasis = true
%!	prm = []
%!	
%!	sizeX = 20
%!	sizeF = sizeX;
%!	numPts = sizeX+1
%!	
%!	vecX0 = randn(sizeX,1) .* exp(expVarCoeff*abs(randn(sizeX,1)));
%!	f0 = max([ exp(expVarCoeff*randn()), 0.0 ]);
%!	matSF = diag(exp(expVarCoeff*randn(sizeX,1)));
%!	matSX = diag(exp(expVarCoeff*randn(sizeX,1)));
%!	matA = randn(sizeX,sizeX) .* exp(expVarCoeff*randn(sizeX,sizeX));
%!	matB = matSF * matA * matSX;
%!	matH = matA' * matA;
%!	condH = cond(matH)
%!	
%!	funchD = @(x)( x - vecX0 );
%!	funchHDOfD = @(d)( matH*d );
%!	funchFSmoothOfD = @(d)( f0 + 0.5*sum( d .* funchHDOfD(d), 1 ) );
%!	funchGSmoothOfD = @(d)( funchHDOfD(d) );
%!	funchFSmooth = @(x)( funchFSmoothOfD(funchD(x)) );
%!	funchGSmooth = @(x)( funchGSmoothOfD(funchD(x)) );
%!	funchF = @(x)( funchFSmooth(x) .* (1.0 + fNoiseLevel*randn([1,size(x,2)])) );
%!	funchG = @(x)( funchGSmooth(x) .* (1.0 + gNoiseLevel*randn(size(x))) );
%!	
%!	matX = randn( sizeX, numPts ) .* exp(expVarCoeff*abs(randn(sizeX,numPts)));
%!	if (forceOrthogonalBasis)
%!		matX(:,1:numPts) = full(eye(sizeX,numPts));
%!	endif
%!	rvecF = funchF(matX);
%!	matG = funchG(matX);
%!	condX = cond(matX)
%!	
%!	[ f_fit, vecG_fit, matH_fit, datOut ] = hessfit( matX, rvecF, matG, prm );
%!	vecXNewton = -(matH_fit\vecG_fit);
%!	
%!	vecG0_true = funchGSmooth(zeros(sizeX,1));
%!	relresG = norm(vecG_fit-vecG0_true) / sqrt( sumsq(vecG_fit)  + sumsq(vecG0_true) )
%!	relresH = sqrt( sum(sum((matH_fit-matH).^2)) ) / sqrt( sum(sum(matH_fit.^2)) + sum(sum(matH.^2)) )
%!	
%!	dist0 = norm(vecX0);
%!	minDistAmongPts = sqrt(min(sum( (matX - vecX0).^2, 1 )));
%!	distNewton = norm( vecXNewton - vecX0 );
%!	dist = [ dist0, minDistAmongPts, distNewton ]
%!	
%!	toc()
