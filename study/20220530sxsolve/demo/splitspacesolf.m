function [ vecX, datOut ] = splitspacesolf( funchMatAProd, vecB, sizeX, prm=[] )
	matVR = mygetfield( prm, "matVR", [] );
	matWR = mygetfield( prm, "matWR", [] );
	tol = mygetfield( prm, "tol", 1e-4 );
	matVL = [];
	matWL = [];
	%msg( __FILE__, __LINE__, "----------" );
	%
	while (1)
		sizeR = size(matVR,2);
		sizeL = size(matVL,2);
		assert( sizeR <= sizeX );
		assert( sizeL <= sizeX );
		if ( 0 == sizeR )
			vecYR = [];
			vecRhoR = vecB;
			resR = norm(vecRhoR);
		else
			vecYR = (matWR'*matWR)\(matWR'*vecB);
			vecRhoR = vecB - matWR*vecYR;
			resR = norm(vecRhoR);
		endif
		%
		if ( 0 == sizeL )
			vecYL = [];
			vecRhoL = vecB;
			resL = norm(vecRhoL);
		else
			vecYL = (matWL'*matWL)\(matWL'*vecB);
			vecRhoL = vecB - matWL*vecYL;
			resL = norm(vecRhoL);
		endif
		%res = [ resR, resL ]
		%
		if ( resL <= tol*norm(vecB) )
			break;
		elseif ( sizeL == sizeX )
			msg( __FILE__, __LINE__, "ALGORITHM BREAKDOWN: Reached full space without convergence." );
			break;
		endif
		%
		%if ( resR > 0.1 * resL && mod(sizeL,3) ~= 0 )
		if ( ( resR > 0.1 * resL && mod(sizeL,3) ~= 0 ) ...
		  || ( resR > 0.9 * resL ) )
			vecU = vecRhoL; % Would apply extra precon here.
			vecV = __orth( vecU, matVR );
			if ( norm(vecV) >= 0.5 )
				%msg( __FILE__, __LINE__, "Expanding." );
				vecW = funchMatAProd( vecV );
				matVR = [ matVR, vecV ];
				matWR = [ matWR, vecW ];
				matVL = [ matVL, vecV ];
				matWL = [ matWL, vecW ];
				continue;
			endif
			%msg( __FILE__, __LINE__, "Cannot expand." );
		endif
		%
		vecU = matVR*vecYR;
		vecV = __orth( vecU, matVL );
		vecV = matVR*(matVR'*vecV); % Help with numerical stability.
		if ( norm(vecV) >= 0.5 )
			vecV /= norm(vecV);
			%
			%msg( __FILE__, __LINE__, "Pulling." );
			vecW = funchMatAProd( vecV );
			vecY = matVR'*vecV;
			assert( norm(vecY)>0.0 );
			vecYHat = vecY/norm(vecY);
			matWR = matWR + ( vecW - matWR*vecYHat )*(vecYHat');
			matVL = [ matVL, vecV ];
			matWL = [ matWL, vecW ];
			continue;
		endif
		%
		vecU = vecRhoL; % Would apply extra precon here.
		vecV = __orth( vecU, matVR );
		if ( norm(vecV) >= 0.5 )
			%msg( __FILE__, __LINE__, "Expanding." );
			vecW = funchMatAProd( vecV );
			matVR = [ matVR, vecV ];
			matWR = [ matWR, vecW ];
			matVL = [ matVL, vecV ];
			matWL = [ matWL, vecW ];
			continue;
		endif
		%
		msg( __FILE__, __LINE__, "Could not find a valid action." );
		break;
	endwhile
	%
	vecX = matVL*vecYL;
	datOut.fevalCount = size(matVL,2);
	datOut.vecY = vecYL;
	datOut.matVL = matVL;
	datOut.matWLs = matWL;
	datOut.matVR = matVR;
	datOut.matWR = matWR;
return;
endfunction



function vecV = __orth( vecU, matV, tol=sqrt(eps) )
	u0 = norm(vecU);
	if (0.0==u0)
		vecV = zeros(size(vecU));
		return;
	elseif (isempty(matV))
		vecV = vecU/u0;
		return;
	elseif ( size(matV,2) >= size(matV,1) )
		vecV = zeros(size(vecU));
		return;
	endif
	vecV = vecU;
	for n=1:2
		vecV -= matV*(matV'*vecV);
		v = norm(vecV);
		if ( v <= tol*u0 )
			vecV(:) = 0.0;
			return;
		else
			vecV /= v;
		endif
	endfor
	return;
endfunction
