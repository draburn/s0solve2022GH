function omegaE = elliptical_omegaBall_evalOmega( funchOmega, vecXC, matS, vecX )
	vecR = vecX - vecXC;
	normR = norm(vecR);
	if ( 0.0 == normR )
		omegaE = funchOmega( vecX );
		return;
	end
	omegaE = funchOmega( vecXC + (matS*vecR/normR) );
return
