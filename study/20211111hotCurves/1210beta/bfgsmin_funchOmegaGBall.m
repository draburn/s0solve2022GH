function [ omega, vecG ] = bfgsmin_funchOmegaGBall( vecX, vecXC, bigR, matS, funchOmega, funchG )
	vecD = vecX - vecXC;
	if ( isempty(matS) )
		s = norm(vecD);
	else
		s = norm(matS*vecD);
	end
	if ( s <= bigR )
		omega = funchOmega( vecX );
		vecG = funchG( vecX );
		return;
	end
	%
	pullCoeff = 0.01;
	vecXSurf = vecXC + (vecD*bigR/s);
	omega = funchOmega( vecXSurf ) + 0.5*pullCoeff*(norm(vecX - vecXSurf)^2);
	vecG0 = funchG( vecXSurf );
	if (isempty(matS))
		vecT = vecD/(s^2);
	else
		vecT = matS'*(matS*vecD)/(s^2);
	end
	vecG = ( vecG0 - vecT*(vecD'*vecG0) )*(bigR/s) + pullCoeff*( vecX - vecXSurf );
	return;
endfunction
