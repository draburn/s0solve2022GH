% Function...
%  solve for vecX: matA * vecX = vecB.
%  More important than vecX, however, is matV and matW;
%   matV provides an orthonormal basis,
%   and matW = matA*matV;
%  This code is simple but slow.

function [ vecX, datOut ] = linsolf_update( funchMatAProd, vecB, prm = [] )
	time0 = time();
	fevalCount = 0;
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	fevalCount = 0;
	%
	sizeF = size(vecB,1);
	sizeX = sizeF;
	assert( isrealarray(vecB,[sizeF,1]) );
	assert( sizeX == sizeF );
	%
	if ( 0.0 == norm(vecB) )
		vecX = zeros(sizeX,1);
		datOut.vecY = [];
		datOut.matV = zeros(sizeX,0);
		datOut.matW = zeros(sizeF,0);
		datOut.vecZ = zeros(sizeF,0);
		datOut.vecRho = zeros(sizeF,0);
		return;
	endif
	%
	matA = mygetfield(prm,"matA",[]);
	if ( ~isempty(matA) )
		assert( isrealarray(matA,[sizeF,sizeX]) );
	endif
	%
	%
	%
	% Everything past here is iterated upon.
	vecU = vecB;
	if ( ~isempty(matA) )
		vecU = matA\vecU;
	endif
	vecV = __orth( vecU, [], prm );
	assert( isrealarray(vecV,[sizeX,1]) );
	vecW = funchMatAProd( vecV ); fevalCount++;
	assert( isrealarray(vecW,[sizeF,1]) );
	matV = vecV;
	matW = vecW;
	while (1)
		sizeV = size(matV,2);
		assert( reldiff( matV'*matV, eye(sizeV,sizeV) ) < sqrt(eps) );
		%
		vecY = __getY( matW, vecB );
		rho = norm( vecB - matW*vecY ) / norm(vecB);
		tol = mygetfield( prm, "tol", 0.1 );
		if ( rho < tol )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "STRONG SUCCESS: rho < tol." );
			break;
		endif
		iterMax = mygetfield( prm, "iterMax", sizeX );
		if ( sizeV >= iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: sizeV >= iterMax." );
			break;
		endif
		%
		vecU = vecW;
		if ( ~isempty(matA) )
			vecU = vecB - matW*vecY; % Isn't this dangerous, because we could get the same vector repeatedly?
			matA += (vecW-matA*vecV)*(vecV'); % THIS IS RATHER HARMFUL?
			vecU = matA\vecU;
		endif
		vecV = __orth( vecU, matV, prm );
		if ( isempty(vecV) )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: Data became linearly dependent." );
			break;
		endif
		assert( isrealarray(vecV,[sizeX,1]) );
		vecW = funchMatAProd( vecV ); fevalCount++;
		assert( isrealarray(vecW,[sizeF,1]) );
		matV = [ matV, vecV ];
		matW = [ matW, vecW ];
	endwhile
	vecX = matV*vecY;
	if ( 2 <= nargout )
		datOut.fevalCount = fevalCount;
		datOut.vecY = vecY;
		datOut.matV = matV;
		datOut.matW = matW;
		datOut.vecZ = matW*vecY;
		datOut.vecRho = vecB - datOut.vecZ;
		datOut.matA = matA;
	endif
return;
end


function [ vecV ] = __orth( vecU, matVPrev, prm=[] )
	uSq = sumsq(vecU);
	if ( uSq == 0.0 )
		vecV = [];
		return;
	endif
	%
	%
	vecV = vecU/sqrt(uSq);
	if (isempty(matVPrev))
		return;
	endif
	%
	numPasses = mygetfield( prm, "numPasses", 2 );
	for n=1:numPasses
		vecV -= matVPrev*(matVPrev'*vecV);
		vSq = sumsq(vecV);
		relTolVSq = mygetfield( prm, "relTolVSq", eps );
		if ( vSq <= relTolVSq*uSq )
			vecV = [];
			return;
		endif
		vecV /= sqrt(vSq);
	endfor
return;
endfunction


function vecY = __getY( matW, vecB, prm=[] )
	matH = matW'*matW;
	sizeV = size(matW,2);
	hNorm = sqrt(sum(sumsq(matH))/sizeV);
	assert( 0 ~= hNorm );
	cholTol = mygetfield( prm, "cholTol", 1e-6 );
	[ matR, cholFlag ] = chol(matH);
	if ( 0==cholFlag && min(diag(matR)) > cholTol*max(abs(diag(matR))) )
		matHRegu = matH;
	else
		matHRegu = matH + hNorm*sqrt(eps)*matIX;
		[ matR, cholFlag ] = chol(matHRegu);
		assert( 0 == cholFlag );
		assert( min(diag(matR)) > max(abs(diag(matR)))*cholTol );
	endif
	vecY = matR \ ( matR' \ (matW'*vecB) );
return;
endfunction
