function [ vecX, datOut ] = semisplitspacelinsolf( funchMatAProd, vecB, vecX0, prm=[] )
	matVR = mygetfield( prm, "matVR", [] );
	matWR = mygetfield( prm, "matWR", [] );
	tol = mygetfield( prm, "tol", 1e-4 );
	sizeX = size(vecX0,1);
	matVL = [];
	matWL = [];
	doproglog = false;
	msgif( doproglog, __FILE__, __LINE__, "----------" );
	sizeF = sizeX;
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
			if ( reldiff( matVR'*matVR, eye(sizeR,sizeR) ) > sqrt(eps) )
				%matVR'*matVR
				%matWR'*matWR
				size(matVR)
				assert( reldiff( matVR'*matVR, eye(sizeR,sizeR) ) < sqrt(eps) );
			endif
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
			if (0)
			if ( reldiff( matVL'*matVL, eye(sizeL,sizeL) ) > sqrt(eps) )
				matVL'*matVL
				matWL'*matWL
				assert( reldiff( matVL'*matVL, eye(sizeL,sizeL) ) < sqrt(eps) );
			endif
			if ( rcond(matWL'*matWL) < eps )
				msg( __FILE__, __LINE__, sprintf("rcond = %g.",rcond(matWL'*matWL)) );
				matVL'*matVL
				matWL'*matWL
				error("WL");
			endif
			endif
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
		%
		%if ( resR > 0.1 * resL && mod(sizeL,3) ~= 0 )
		%if ( ( resR > 0.1 * resL && mod(sizeL,3) ~= 0 ) ...
		%  || ( resR > 0.9 * resL ) )
		c = [ 0.01, 0.1, 0.5, 0.9, 0.99 ](1+mod(sizeL,5));
		if ( resR > c*resL )
			vecU = __applyPrecon( vecRhoL, prm );
			%%%vecV = __orth( vecU, matVL );
			vecV = __orth( vecU, matVR ); % VIOLATES "SEMI-SPLIT-SPACE" CONCEPT, BUT, MAY BE A GOOD IDEA?
			if ( norm(vecV) >= 0.5 )
				msgif( doproglog, __FILE__, __LINE__, "Expanding." );
				vecW = funchMatAProd( vecV );
				if ( sizeR == sizeX )
					vecY = matVR'*vecV;
					matWR += ( vecW - matWR*vecY )*(vecY');
				elseif ( norm(matVR'*vecV) > 100.0*eps )
					vecVPerp = __orth(vecV,matVR);
					assert( norm(vecVPerp) > 0.5 );
					matVR = [ matVR, vecVPerp ];
					vecY = matVR'*vecV;
					matWR = [ matWR, zeros(sizeF,1) ];
					matWR += ( vecW - matWR*vecY) * (vecY');
				else
					matVR = [ matVR, vecV ];
					matWR = [ matWR, vecW ];
				endif
				matVL = [ matVL, vecV ];
				matWL = [ matWL, vecW ];
				continue;
			endif
			msgif( doproglog, __FILE__, __LINE__, "Cannot expand." );
		endif
		%
		vecU = matVR*vecYR;
		vecV = __orth( vecU, matVL );
		vecV = matVR*(matVR'*vecV); % Help with numerical stability.
		if ( norm(vecV) >= 0.5 )
			vecV /= norm(vecV);
			%
			msgif( doproglog, __FILE__, __LINE__, "Pulling." );
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
		vecU = __applyPrecon( vecRhoL, prm );
		%%%vecV = __orth( vecU, matVL );
		vecV = __orth( vecU, matVR ); % VIOLATES "SEMI-SPLIT-SPACE" CONCEPT, BUT, MAY BE A GOOD IDEA?
		if ( norm(vecV) >= 0.5 )
			msgif( doproglog, __FILE__, __LINE__, "Expanding (2nd pass)." );
			vecW = funchMatAProd( vecV );
			if ( sizeR == sizeX )
				vecY = matVR'*vecV;
				matWR += ( vecW - matWR*vecY )*(vecY');
			elseif ( norm(matVR'*vecV) > 100.0*eps )
				vecVPerp = __orth(vecV,matVR);
				assert( norm(vecVPerp) > 0.5 );
				matVR = [ matVR, vecVPerp ];
				vecY = matVR'*vecV;
				matWR = [ matWR, zeros(sizeF,1) ];
				matWR += ( vecW - matWR*vecY) * (vecY');
			else
				matVR = [ matVR, vecV ];
				matWR = [ matWR, vecW ];
			endif
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
	datOut.matWL = matWL;
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

function vecU = __applyPrecon( vecRho, prm )
	%vecU = vecRho;
	%
	%s = sqrt( sum(sumsq(matWR)) / sum(sumsq(matVR)) );
	%sizeX = size(matVR,1);
	%vecU = ( s*eye(sizeX,sizeX) + ( matWR - matVR ) * (matVR') ) \ vecRho;
	vecU = prm.matP*vecRho;
	return;
endfunction
