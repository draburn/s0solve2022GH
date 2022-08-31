% 800 simplified: no AP, no coasting, and (perhaps?) no BT.
%  also: match epsFD
% Watch linsolf params.

function [ vecXF, vecFF, datOut ] = findZero_800simple( vecX0, funchF, prm=[] )
	time0 = time();
	fevalCount = 0;
	mydefs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	vecF0 = funchF( vecX0 ); fevalCount++;
	sizeF = size(vecF0,1);
	assert( isrealarray(vecF0,[sizeF,1]) );
	assert( 0~=norm(vecF0) );
	%
	epsFD = mygetfield( prm, "epsFD", 1.0e-4 ); % Was eps^0.3 or eps^0.4.
	fTol = mygetfield( prm, "fTol", sqrt(sizeF)*100.0*eps );
	iterMax = mygetfield( prm, "iterMax", 50 );
	dumpForDebug = false;
	%
	%
	%
	matIX = eye(sizeX,sizeX);
	%
	%
	% Everything past here is iterated on.
	iterCount = 0;
	vecX_best = vecX0;
	vecF_best = vecF0;
	datOut.fNormVals(iterCount+1) = norm(vecF_best);
	datOut.fevalCountVals(iterCount+1) = fevalCount;
	datOut.iterCountVals(iterCount+1) = iterCount;
	%
	vecX = vecX0;
	vecF = vecF0;
	%
	msgif( prm.verbLev >= VERBLEV__INFO, __FILE__, __LINE__, sprintf( ...
	  "  %8.3e, %3d, %3d: %4d / %dx%d;  %8.2e // %8.2e // %8.2e / %8.2e / %8.2e.", ...
	  time()-time0, iterCount, fevalCount, ...
	  0, sizeF, sizeX, ...
	  (fTol^2)/2.0, sumsq(vecF)/2.0, sumsq(vecF)/2.0, sumsq(vecF)/2.0, sumsq(vecF0)/2.0 ) );
	%
	%
	%xf = [ vecX, vecF ]
	funchMatJProd = @(v)( ( funchF(vecX+epsFD*v) - vecF ) / epsFD );
	linsolf_prm = [];
	linsolf_prm.tol = mygetfield( prm, "linsolf_tol", 0.1*sqrt(norm(vecF)/norm(vecF0)) );
	linsolf_prm.tol = mygetfield( prm, "linsolf_tol0", linsolf_prm.tol );
	linsolf_prm = mygetfield( prm, "linsolf_prm", linsolf_prm );
	[ vecSSDeltaN, linsolf_datOut ] = linsolf_simple( funchMatJProd, -vecF, linsolf_prm );
	fevalCount += linsolf_datOut.fevalCount;
	sizeV = size(linsolf_datOut.matV,2);
	matV = linsolf_datOut.matV;
	matW = linsolf_datOut.matW;
	if (dumpForDebug)
	epsFD
	for k=1:sizeV
		k
		vecV = matV(:,k);
		vecW = matW(:,k);
		vwxf = [ vecV, vecW, vecX, vecF, vecX+epsFD*vecV, funchF(vecX+epsFD*vecV) ]
	endfor
	endif
	%
	vecX_next = vecX - matV*( (matW'*matW) \ (matW'*vecF) );
	vecFModel_next = vecF - matW*( (matW'*matW) \ (matW'*vecF) );
	vecF_next = funchF(vecX_next); fevalCount++;
	%
	%
	%
	while (1)
		iterCount++;
		if ( norm(vecF_next) < norm(vecF_best) )
			vecX_best = vecX_next;
			vecF_best = vecF_next;
		endif
		datOut.fNormVals(iterCount+1) = norm(vecF_best);
		datOut.fevalCountVals(iterCount+1) = fevalCount;
		datOut.iterCountVals(iterCount+1) = iterCount;
		%
		msgif( prm.verbLev >= VERBLEV__INFO, __FILE__, __LINE__, sprintf( ...
		  "  %8.3e, %3d, %3d: %4d / %dx%d;  %8.2e // %8.2e // %8.2e / %8.2e / %8.2e.", ...
		  time()-time0, iterCount, fevalCount, ...
		  sizeV, sizeF, sizeX, ...
		  (fTol^2)/2.0, sumsq(vecFModel_next)/2.0, sumsq(vecF_next)/2.0, sumsq(vecF)/2.0, sumsq(vecF0)/2.0 ) );
		%
		%
		%
		if ( norm(vecF_best) <= fTol )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "STRONG SUCCESS: norm(vecF_best) <= fTol." );
			break;
		elseif ( norm(vecF_next) + fTol >= norm(vecF) )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: norm(vecF_next) + fTol >= norm(vecF)." );
			break;
		elseif ( iterCount >= iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterMax." );
			break;
		endif
		vecX = vecX_next;
		vecF = vecF_next;
		%
		%
		funchMatJProd = @(v)( ( funchF(vecX+epsFD*v) - vecF ) / epsFD ); % vecX and vecF have been updated!
		linsolf_prm = [];
		linsolf_prm.tol = mygetfield( prm, "linsolf_tol", 0.1*sqrt(norm(vecF)/norm(vecF0)) );
		linsolf_prm = mygetfield( prm, "linsolf_prm", linsolf_prm );
		[ vecSSDeltaN, linsolf_datOut ] = linsolf_simple( funchMatJProd, -vecF, linsolf_prm );
		fevalCount += linsolf_datOut.fevalCount;
		sizeV = size(linsolf_datOut.matV,2);
		matV = linsolf_datOut.matV;
		matW = linsolf_datOut.matW;
		if (dumpForDebug)
		epsFD
		for k=1:sizeV
			k
			vecV = matV(:,k);
			vecW = matW(:,k);
			vwxf = [ vecV, vecW, vecX, vecF, vecX+epsFD*vecV, funchF(vecX+epsFD*vecV) ]
		endfor
		endif
		%
		%
		vecX_next = vecX - matV*( (matW'*matW) \ ( matW'*vecF ) );
		vecFModel_next = vecF - matW*( (matW'*matW) \ (matW'*vecF) );
		vecF_next = funchF(vecX_next); fevalCount++;
		%xfTrial = [ vecX_next, vecF_next ]
	endwhile
	vecXF = vecX_best;
	vecFF = vecF_best;
	datOut.fevalCount = fevalCount;
	datOut.iterCount = iterCount;
	%
return;
endfunction
