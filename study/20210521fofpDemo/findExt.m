thisFile = "findExt";

% Calculate our "super points".
% These are, in some sense (f')/(f'').
numSuperPts = numPts-2;
for m=1:numSuperPts
	vecX = xVals(m:m+2)';
	matX = [ ones(3,1), vecX, vecX.^2 ];
	vecF = fVals(m:m+2)';
	vecC = matX \ vecF;
	c0 = vecC(1);
	c1 = vecC(2);
	c2 = vecC(3);
	%
	% What should super_xVal be???
	super_xVals(m) = 0.5*( vecX(1) + vecX(3) );
	super_xVals_probablyBad(m) = vecX(2);
	%
	super_hVals(m) = 0.5*c1/c2 + super_xVals(m);
end
%
% From our "super points", find a best-fit bigX and p.
super_matH = [ ones(numSuperPts,1), super_hVals' ];
super_vecRes = [ super_hVals' + super_xVals' ];
super_vecC = super_matH \ super_vecRes;
bigX_initial = super_vecC(1);

if (0)
	bigX_initial = 0.2
	msg( thisFile, __LINE__, sprintf("HACK: Set bigX_initial = %f.", bigX_initial) );
end

bigP_initial = super_vecC(2);
matK_initial = [ ones(numPts,1), abs(xVals'-bigX_initial).^bigP_initial ];
vecF = fVals';
vecBigF_initial = matK_initial \ vecF;
bigF0_initial = vecBigF_initial(1);
bigF1_initial = vecBigF_initial(2);

msg( thisFile, __LINE__, "TODO: Consider that our bigX may be in wrong interval!" );
msg( thisFile, __LINE__, "TODO: Prevent P from going negative?" );
%
%
doLoop = true;
loopCount = 0;
btCount = 0;
vecDelta0 = zeros(2,1);
vecDelta = zeros(2,1);
bigX = bigX_initial;
bigP = bigP_initial;
msg( thisFile, __LINE__, sprintf( ...
  " %3d, %2d,   %10.3e, %10.3e, %10.3e, %10.3e,   %10.3e, %10.3e, %10.3e, %10.3e", ...
  loopCount, btCount, ...
  vecDelta0(1), vecDelta(1), bigX, bigX-secret_bigX, ...
  vecDelta0(2), vecDelta(2), bigP, bigP-secret_bigP ) );
while (doLoop)
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
	doBTLoop = true;
	btCount = 0;
	while (doBTLoop)
		bigX_temp = bigX + vecDelta(1);
		bigP_temp = bigP + vecDelta(2);
		rhoVals_temp = findExt__rho( xVals, fVals, bigX_temp, bigP_temp );
		if ( sum(rhoVals_temp.^2) < sum(rhoVals_0.^2) )
			doBTLoop = false;
		else
			btCount++;
			vecDelta /= 5.0;
			if ( btCount >= 10 )
				doBTLoop = false;
			end
		end
	end
	%
	bigX = bigX + vecDelta(1);
	bigP = bigP + vecDelta(2);
	%
	loopCount++;
	if ( loopCount>=50 )
		doLoop = false;
	end
	msg( thisFile, __LINE__, sprintf( ...
	  " %3d, %2d,   %10.3e, %10.3e, %10.3e, %10.3e,   %10.3e, %10.3e, %10.3e, %10.3e", ...
	  loopCount, btCount, ...
	  vecDelta0(1), vecDelta(1), bigX, bigX-secret_bigX, ...
	  vecDelta0(2), vecDelta(2), bigP, bigP-secret_bigP ) );
end

% Finally: Find corresponding F0 and F1
matK = [ ones(numPts,1), abs(xVals'-bigX).^bigP ];
vecF = fVals';
vecBigF = matK \ vecF;
bigF0 = vecBigF(1);
bigF1 = vecBigF(2);
return;
