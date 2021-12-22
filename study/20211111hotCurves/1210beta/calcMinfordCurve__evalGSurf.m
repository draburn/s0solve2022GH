function vecGSurf = calcMinfordCurve__evalGSurf( funchG, vecXC, bigR, vecX, matS=[], prm=[] )
	%thisFile = "calcMinfordCurve__evalGSurf";
	%
	doChecks = mygetfield(prm,"doChecks",true);
	if ( doChecks )
		sizeX = size(vecXC,1);
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
		vecGSurf = funchG( vecXC );
		return;
	end
	%
	if (isempty(matS))
		vecDHat = vecD/normD;
		vecXSurf = vecXC + (bigR*vecDHat);
		vecG0 = funchG( vecXSurf );
		vecGSurf0 = vecG0 - (vecDHat*(vecDHat'*vecG0));
		vecGSurf = vecGSurf0*(bigR/normD);
	elseif (0)
		% Old/bad version.
		vecDHat = vecD/normD;
		vecGSurf0 = matS'*(funchG( vecXC + (bigR*matS*vecDHat) )*(bigR/normD));
		vecGSurf = vecGSurf0 - (vecDHat*(vecDHat'*vecGSurf0));
	else
		vecSD = matS*vecD;
		normSD = norm(vecSD);
		vecXSurf = vecXC + (bigR*vecD/normSD);
		vecG0 = funchG( vecXSurf );
		vecGSurf0 = vecG0 - (matS'*matS*vecD)*(vecD'*vecG0)/(normSD^2);
		vecGSurf = vecGSurf0*(bigR/normSD);
	end
return;
end
