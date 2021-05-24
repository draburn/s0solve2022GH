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
	%super_xVals(m) = vecX(2);
	%
	super_hVals(m) = 0.5*c1/c2 + super_xVals(m);
end
%
% From our "super points", find a best-fit bigX and p.
super_matH = [ ones(numSuperPts,1), super_hVals' ];
super_vecRes = [ super_hVals' + super_xVals' ];
super_vecC = super_matH \ super_vecRes;
bigX = super_vecC(1);
bigP = super_vecC(2);
%
%
doLoop = true
loopCount = 0;
while (doLoop)
	temp_bigX = bigX
	temp_bigP = bigP
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
	vecDelta = ...
	 [ sum(dRho_dx.^2),       sum(dRho_dx.*dRho_dp); ...
	   sum(dRho_dx.*dRho_dp), sum(dRho_dp.^2) ...
	 ] \ [ sum(dRho_dx.*rhoVals_0); sum(dRho_dp.*rhoVals_0) ];
	%
	% Simple geom backtrack to prevent p going negative.
	while (  bigP + vecDelta(2) < 0 )
		vecDelta *= 0.1/vecDelta(2);
	end
	%
	bigX = bigX + vecDelta(1)
	bigP = bigP + vecDelta(2)
	%
	loopCount++;
	echo__loopCount = loopCount
	if (loopCount>=20)
		doLoop = false;
	end
end

% Finally: Find corresponding F0 and F1
matK = [ ones(numPts,1), abs(xVals'-bigX).^bigP ];
vecF = fVals';
vecBigF = matK \ vecF;
bigF0 = vecBigF(1);
bigF1 = vecBigF(2);
return;


%
% Next: Nonlinear / iterative-linear solve for bigX and p.
"TO DO!"
bigX = bigX
bigP = bigP
%
% Finally: Find corresponding F0 and F1
matK = [ ones(numPts,1), abs(xVals'-bigX).^bigP ];
vecF = fVals';
vecBigF = matK \ vecF;
bigF0 = vecBigF(1);
bigF1 = vecBigF(2);
%
"Consider doing a 4-param nonlinear / iterative-linear optimization."
bigX = bigX
bigP = bigP
bigF0 = bigF0
bigF1 = bigF1