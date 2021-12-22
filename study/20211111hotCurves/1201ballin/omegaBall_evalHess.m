function matHBall = omegaBall_evalHess( funchG, funchH, vecXC, bigR, bigC, vecX )
	sizeX = size(vecXC,1);
	assert( isrealarray(vecXC,[sizeX,1]) );
	assert( isrealscalar(bigR) );
	assert( 0.0 < bigR );
	assert( isrealscalar(bigC) );
	assert( isrealarray(vecX,[sizeX,1]) );
	%
	vecR = vecX - vecXC;
	normR = norm(vecR);
	if ( normR < bigR )
		matHBall = funchH(vecX);
		return;
	end
	%
	vecRHat = vecR/normR;
	vecG = funchG( vecXC + (vecR*bigR/normR) );
	matH = funchH( vecXC + (vecR*bigR/normR) );
	matI = eye(sizeX,sizeX);
	%
	% I hate my notation. =(
	%matHBall = bigR/(normR^2) * (  ...
	%   3.0*(vecRHat'*vecG)*vecRHat*(vecRHat') ...
	%   - (vecRHat'*vecG)*matI ...
	%   - vecRHat*(vecG') ...
	%   - vecG*(vecRHat')  ) ...
	% + ((bigR/normR)^2)*( matI - vecRHat*(vecRHat'))*matH*( matI - vecRHat*(vecRHat')) ...
	% + 4.0*bigC*( (normR^2-bigR^2)*matI + 2.0*vecR*(vecR') );
	%
	mat1a = 3.0*(vecRHat'*vecG)*vecRHat*(vecRHat') - (vecRHat'*vecG)*matI - vecRHat*(vecG') - vecG*(vecRHat');
	mat1b = mat1a*bigR/(normR^2);
	mat2a = ( matI - vecRHat*(vecRHat') ) * matH * ( matI - vecRHat*(vecRHat') );
	mat2b = mat2a*(bigR^2)/(normR^2);
	mat3 = 4.0*bigC*( (normR^2-bigR^2)*matI + 2.0*vecR*(vecR') );
	matHBall = mat1b + mat2b + mat3;
return
