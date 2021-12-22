function matHE = elliptical_omegaBall_evalHess( funchGrad, funchHess, vecXC, matS, vecX )
	sizeX = size(vecXC,1);
	assert( isrealarray(vecXC,[sizeX,1]) );
	assert( isrealarray(matS,[sizeX,sizeX]) );
	assert( isrealarray(vecX,[sizeX,1]) );
	%
	vecR = vecX - vecXC;
	normR = norm(vecR);
	if ( 0.0 == normR )
		matHE = funchH(vecX);
		return;
	end
	%
	vecRHat = vecR/normR;
	vecSTG = matS'*funchGrad( vecXC + (matS*vecR/normR) );
	matSTHS = matS'*funchHess( vecXC + (matS*vecR/normR) )*matS;
	matI = eye(sizeX,sizeX);
	%
	mat1 = 3.0*(vecRHat'*vecSTG)*vecRHat*(vecRHat') - (vecRHat'*vecSTG)*matI - vecRHat*(vecSTG') - vecSTG*(vecRHat');
	mat2 = ( matI - vecRHat*(vecRHat') ) * matSTHS * ( matI - vecRHat*(vecRHat') );
	matHE = ( mat1 + mat2 ) / (normR^2);
return
