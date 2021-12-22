function vecGamma = findSphereExt__calcGamma( funchG, vecXC, bigR, vecX, prm=[] )
	funchBeta = mygetfield( prm, "funchBeta", @(rSq)(eps*( rSq - (bigR^2) )) );
	vecR = vecX - vecXC;
	normRSq = sum(vecR.^2);
	assert( 0.0 ~= normRSq );
	vecG = funchG(vecXC + (bigR*vecR/sqrt(normRSq)));
	vecGamma = vecG - ((vecR*(vecR'*vecG))/normRSq) + (vecR*funchBeta(normRSq));
return;
end
