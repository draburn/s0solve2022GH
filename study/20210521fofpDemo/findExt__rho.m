function rhoVals = findExt__f( xVals, fVals, bigX, bigP )
	thetaVals = abs( xVals - bigX ).^bigP;
	numPts = max(size(xVals));
	sumThetaVals = sum(thetaVals);
	vecF = ....
	  [ numPts, sumThetaVals; sumThetaVals, sum(thetaVals.^2) ] ...
	  \ [ sum(fVals); sum(thetaVals.*fVals) ];
	bigF0 = vecF(1);
	bigF1 = vecF(2);
	rhoVals = bigF0 + bigF1 * abs( xVals - bigX ).^bigP;
return;
end
