thisFile = "findExt__loop";
%
loopCount = 0;
btCount = 0;
vecDelta0 = zeros(2,1);
vecDelta = zeros(2,1);
%
msg( thisFile, __LINE__, sprintf( ...
  "  %3d, %2d, %9.3e,   %10.3e, %10.3e, %10.3e, %10.3e,   %10.3e, %10.3e, %10.3e, %10.3e", ...
  loopCount, btCount, res, ...
  vecDelta0(1), vecDelta(1), bigX, bigX-secret_bigX, ...
  vecDelta0(2), vecDelta(2), bigP, bigP-secret_bigP ) );
while (1)
	%
	epsBigX = sqrt(eps)*abs(max(xVals)-min(xVals));
	epsBigP = sqrt(eps);
	rhoVals_0  = findExt__rho( xVals, fVals, bigX,         bigP         );
	rhoVals_px = findExt__rho( xVals, fVals, bigX+epsBigX, bigP         );
	rhoVals_mx = findExt__rho( xVals, fVals, bigX-epsBigX, bigP         );
	rhoVals_pp = findExt__rho( xVals, fVals, bigX,         bigP+epsBigP );
	rhoVals_mp = findExt__rho( xVals, fVals, bigX,         bigP-epsBigP );
	%
	dRho_dx = 0.5 * ( rhoVals_px - rhoVals_mx ) / epsBigX;
	dRho_dp = 0.5 * ( rhoVals_pp - rhoVals_mp ) / epsBigP;
	%
	vecDelta0 = -...
	 [ sum(dRho_dx.^2),       sum(dRho_dx.*dRho_dp); ...
	   sum(dRho_dx.*dRho_dp), sum(dRho_dp.^2) ...
	 ] \ [ sum(dRho_dx.*rhoVals_0); sum(dRho_dp.*rhoVals_0) ];
	vecDelta = vecDelta0;
	%
	btCount = 0;
	if ( bigP + vecDelta(2) <= sqrt(eps) )
		vecDelta *= 0.2*bigP/norm(vecDelta);
		btCount++;
	end
	%
	maxBTCount = 10;
	while (1)
		bigX_temp = bigX + vecDelta(1);
		bigP_temp = bigP + vecDelta(2);
		rhoVals_temp = findExt__rho( xVals, fVals, bigX_temp, bigP_temp );
		if ( 0.5*sum(rhoVals_temp.^2) < 0.5*sum(rhoVals_0.^2) )
			break;
		end
		btCount++;
		vecDelta /= 5.0;
		if ( btCount >= maxBTCount )
			break;
		end
	end
	if ( btCount >= maxBTCount )
		msg( thisFile, __LINE__, "Failed to decrease residual." );
		break;
	end
	%
	bigX = bigX + vecDelta(1);
	bigP = bigP + vecDelta(2);
	loopCount++;
	res = 0.5*sum(rhoVals_temp.^2);
	msg( thisFile, __LINE__, sprintf( ...
	  " %3d, %2d, %9.3e,   %10.3e, %10.3e, %10.3e, %10.3e,   %10.3e, %10.3e, %10.3e, %10.3e", ...
	  loopCount, btCount, res, ...
	  vecDelta0(1), vecDelta(1), bigX, bigX-secret_bigX, ...
	  vecDelta0(2), vecDelta(2), bigP, bigP-secret_bigP ) );	
	if ( loopCount >= maxLoopCount )
		break;
	end
	if ( res <= resTarget )
		break;
	end
end
