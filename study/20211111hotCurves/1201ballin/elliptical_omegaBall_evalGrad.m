function vecGradE = elliptical_omegaBall_evalGrad( funchGrad, vecXC, matS, vecX )
	vecR = vecX - vecXC;
	normR = norm(vecR);
	if ( 0.0 == normR )
		vecGradE = funchGrad(vecX);
		return;
	end
	vecRHat = vecR/normR;
	vecSTG = matS'*funchGrad( vecXC + (matS*vecRHat) );
	vecGradE = ( vecSTG - vecRHat*(vecRHat'*vecSTG) ) / normR;
return
