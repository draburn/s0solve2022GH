function [ rhoVals, bigF0, bigF1, omega, errFlag ] = extFit__calcAtPt_turbo( s, p, xVals, fVals, dVals )
	yVals = abs(xVals-s).^p;
	if (~isrealarray(yVals,size(xVals)))
		rhoVals = -1.0;
		bigF0 = 0.0;
		bigF1 = 0.0;
		omega = -1.0;
		errFlag = true;
		return;
		% For example, 47^209 = "Inf".
	end
	matY = [ ones(size(yVals))', yVals' ];
	matD = diag(dVals);
	matA = matD*matY;
	matB = matA'*matA;
	if ( eps >= rcond(matB) )
		rhoVals = -1.0;
		bigF0 = 0.0;
		bigF1 = 0.0;
		omega = -1.0;
		errFlag = true;
		return;
	end
	vecC = matB\(matA'*matD*(fVals'));
	bigF0 = vecC(1);
	bigF1 = vecC(2);
	rhoVals = bigF0 + (bigF1*yVals) - fVals;
	omega = 0.5*sum(( dVals.*rhoVals ).^2);
	errFlag = false;
return;
end
