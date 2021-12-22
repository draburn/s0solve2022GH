function vecGBall = omegaBall_evalGrad( funchG, vecXC, bigR, bigC, vecX )
	vecR = vecX - vecXC;
	normR = norm(vecR);
	if ( normR < bigR )
		vecGBall = funchG(vecX);
		return;
	end
	vecRHat = vecR/normR;
	vecG = funchG( vecXC + (vecR*bigR/normR) );
	vecGBall = (bigR/normR)*( vecG - vecRHat*(vecRHat'*vecG)) ...
	  + 4.0*bigC*(normR^2-bigR^2)*vecR;
return
