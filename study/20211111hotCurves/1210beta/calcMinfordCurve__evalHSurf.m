function matHSurf = calcMinfordCurve__evalHSurf( funchG, funchH, vecXC, bigR, vecX, matS=[], prm=[] )
	%thisFile = "calcMinfordCurve__evalHSurf";
	%
	sizeX = size(vecXC,1);
	doChecks = mygetfield(prm,"doChecks",true);
	if ( doChecks )
		assert( isrealarray(vecXC,[sizeX,1]) );
		assert( isrealscalar(bigR) );
		assert( bigR > 0.0 );
		assert( isrealarray(vecX,[sizeX,1]) );
		if (~isempty(matS))
			assert( isrealarray(matS,[sizeX,sizeX]) );
		end
	end
	%
	vecD = vecX - vecXC;
	normD = norm(vecD);
	if ( 0.0 == normD )
		matHSurf = funchH( vecXC );
		return;
	end
	%
	vecDHat = vecD/normD;
	matIX = eye(sizeX,sizeX);
	if (isempty(matS))
		vecXSurf = vecXC + (bigR*vecDHat);
		vecG = funchG( vecXSurf );
		matH = funchH( vecXSurf );
		mat1 = 3.0*(vecG'*vecDHat)*vecDHat*(vecDHat') - (vecG'*vecDHat)*matIX - vecDHat*(vecG') - vecG*(vecDHat');
		mat2 = ( matIX - vecDHat*(vecDHat') ) * matH * ( matIX - vecDHat*(vecDHat') );
		matHSurf = (mat1*(bigR/(normD^2))) + (mat2*((bigR/normD)^2));
	elseif (0)
		% Old/bad version.
		vecXSurf = vecXC + (bigR*matS*vecDHat);
		vecG = funchG( vecXSurf );
		matH = funchH( vecXSurf );
		mat1 = 3.0*(vecG'*matS*vecDHat)*vecDHat*(vecDHat') - (vecG'*matS*vecDHat)*matIX - vecDHat*(vecG'*matS) - (matS'*vecG)*(vecDHat');
		mat2 = ( matIX - vecDHat*(vecDHat') ) * matS' * matH * matS * ( matIX - vecDHat*(vecDHat') );
		matHSurf = (mat1*(bigR/(normD^2))) + (mat2*((bigR/normD)^2));
	else
		vecSD = matS*vecD;
		normSD = norm(vecSD);
		vecSTSD = matS'*vecSD;
		vecXSurf = vecXC + (bigR*vecD/normSD);
		vecG = funchG( vecXSurf );
		matH = funchH( vecXSurf );
		matA = (bigR/normSD)*( matIX - (vecSTSD*(vecD'))/(normSD^2) );
		mat1 = matA*matH*(matA');
		mat2 = 3.0*bigR*(vecD'*vecG)*(vecSTSD*(vecSTSD'))/(normSD^5);
		mat3 = (bigR/(normSD^3))*( vecSTSD*(vecG') + vecG*(vecSTSD') );
		mat4 = (matS'*matS)*(bigR*(vecD'*vecG)/(normSD^3));
		matHSurf = mat1 + mat2 - mat3 - mat4;
	end
	matHSurf = (matHSurf'+matHSurf)/2.0;
return;
end
