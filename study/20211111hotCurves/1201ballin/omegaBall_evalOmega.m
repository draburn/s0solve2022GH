function omega = omegaBall_evalOmega( funchOmega, vecXC, bigR, bigC, vecX )
	vecR = vecX - vecXC;
	normR = norm(vecR);
	if ( normR < bigR )
		omega = funchOmega( vecX );
	else
		omega = funchOmega( vecXC + (vecR*bigR/normR) ) + (bigC*(normR^2-bigR^2)^2);
	end
return
