function omega = bfgsmin_funchOmegaBall( vecX, vecXC, bigR, matS, funchOmega )
	vecD = vecX - vecXC;
	if ( isempty(matS) )
		s = norm(vecD);
	else
		s = norm(matS*vecD);
	end
	if ( s <= bigR )
		omega = funchOmega( vecX );
		return;
	end
	vecXBall = vecXC + (vecD*bigR/s);
	vecQ = vecX - vecXBall;
	omega = funchOmega( vecXBall ) + 0.1*(vecQ'*vecQ) + 1e-8*(vecD'*vecD);
	return;
endfunction
