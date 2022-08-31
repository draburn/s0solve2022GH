function [ vecX, datOut ] = linsolf_simple( funchMatAProd, vecB, prm = [] )
	mydefs;
	time0 = time();
	verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	datOut = [];
	fevalCount = 0;
	%
	tol = mygetfield( prm, "tol", 0.1 );
	sz = size(vecB,1);
	%
	vecV = __orth( vecB, [], prm );
	if ( 0.0 == norm(vecV) )
		error( "Could not create new subspace vector." );
	endif
	vecW = funchMatAProd( vecV ); fevalCount++;
	%vw = [ vecV, vecW ]
	matV = vecV;
	matW = vecW;
	%
	while (1)
		vecX = ( (matW'*matW) \ (matW'*vecB) );
		vecRho = vecB - matW*vecX;
		if ( norm(vecRho) < tol*norm(vecB) )
			break;
		endif
		if ( size(matV,2) == sz )
			error( "Reached full space without converging." );
		endif
		%%%vecV = __orth( vecW, matV, prm );
		vecV = __orth( vecRho, matV, prm );
		if ( 0.0 == norm(vecV) )
			error( "Could not create new subspace vector." );
		endif
		vecW = funchMatAProd( vecV ); fevalCount++;
		%vw = [ vecV, vecW ]
		matV = [ matV, vecV ];
		matW = [ matW, vecW ];
	endwhile
	%
	datOut.fevalCount = fevalCount;
	datOut.matV = matV;
	datOut.matW = matW;
	return;
endfunction


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
