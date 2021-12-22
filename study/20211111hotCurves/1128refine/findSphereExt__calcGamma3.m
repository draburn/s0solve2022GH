function vecGamma = findSphereExt__calcGamma( funchG, vecXC, bigR, vecX )
	vecR = vecX - vecXC;
	normRSq = sum(vecR.^2);
	assert( 0.0 ~= normRSq );
	normR = sqrt(normRSq);
	vecG = funchG(vecXC + (bigR*vecR/normR));
	vecGamma = (vecG - ((vecR*(vecR'*vecG))/normRSq))*(bigR/normR) ...
	  + 2.0*(normRSq-bigR^2)*vecR;
return;
end
