	clear;
	numFigs = 0;
	setprngstates(0);
	sizeX = 200;
	sizeF = 100;
	%
	vecF = randn(sizeF,1);
	matJ = randn(sizeF,sizeX);
	if (1)
		vecPhi = randn(sizeX,1);
		vecPhi /= norm(vecPhi);
		matJ -= (matJ*vecPhi)*(vecPhi');
	endif
	matB_gen = randn(sizeF,sizeX);
	matB_unscaled = matB_gen'*matB_gen;
	%matB_unscaled = diag(diag(matJ'*matJ));
	bMax_unscaled = sqrt(2.0);
	%
	%
	matH = matJ'*matJ;
	vecG = matJ'*vecF;
	%
	% Regularize things?
	matB_unscaled += sqrt(eps)*max(diag(matB_unscaled))*eye(sizeX,sizeX);
	matH += sqrt(eps)*max(diag(matH))*eye(sizeX,sizeX);
	%
	% Scale B?
	%hobScale =  sqrt(max(diag(matH))) / max(diag(matB_unscaled));
	%hobScale =  1.0e-4*sqrt(max(diag(matH))) / max(diag(matB_unscaled))
	%bPrime0 = norm(matB_unscaled*( (matB_unscaled'*matB_unscaled) \ vecG ))
	%b1 = norm(matB_unscaled*( matH \ vecG ))
	%hobScale =  sqrt(bPrime0/b1)
	%hobScale = sqrt( norm(matH\vecG) / norm(matH\(matB_unscaled'*(matB_unscaled*(matH\vecG)))) )
	vecHUG = matH\vecG;
	vecSHUG = matB_unscaled' * ( matB_unscaled*vecHUG );
	if (0)
	vecHUSHUG = matH \ vecSHUG;
	[ norm( vecSHUG ), norm(vecHUSHUG) ]
	shugthushug = vecSHUG' * vecHUSHUG
	hobScale2 = norm( matB_unscaled*vecHUG ) / sqrt( shugthushug )
	a = sumsq( matB_unscaled * vecHUSHUG )
	b = 2.0*( vecSHUG' * vecHUSHUG )
	c = sumsq( matB_unscaled * vecHUG ) * 3.0/4.0
	discrim = b^2 - 4 * a * c
	a*(hobScale2^4)
	b*(hobScale2^2)
	return;
	hobScale = hobScale2
	else
		matS_foo = matB_unscaled'*matB_unscaled;
		matS_foo *= sqrt(eps) * max(diag(matH)) / max(diag(matS_foo));
		matR1 = chol( matH + matS_foo );		
		matR2 = chol( matH + 2.0*matS_foo );
		vecRTUSHUG = (2.0*( matR1' \ vecSHUG )) - ( matR2' \ vecSHUG );
		%hobScale = norm( matB_unscaled*vecHUG ) / norm( vecRTUSHUG )
		hobScale = 1.0;
		%matR = chol( matH );
		%vecRTUSHUG = matR' \ vecSHUG;
		%hobScale = norm( matB_unscaled*vecHUG ) / norm( vecRTUSHUG )
		%return;
	endif
	matB_scaled = matB_unscaled * hobScale;
	bMax_scaled = bMax_unscaled * hobScale;
	%rcond( matH + matB_scaled'*matB_scaled )
	%
	%
	prm = [];
	prm.bTol = 0.001*bMax_scaled;
	[ vecY, vecYPrime, b, bPrime ] = findLevPt( vecG, matH, bMax_scaled, matB_scaled, prm );
	[ norm(matB_scaled*vecY), bMax_scaled, norm(matB_scaled*vecY) - bMax_scaled ]
	[ norm(matB_unscaled*vecY), bMax_unscaled, norm(matB_unscaled*vecY) - bMax_unscaled ]
