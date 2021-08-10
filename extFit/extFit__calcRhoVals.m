function [ rhoVals, errFlag ] = extFit__calcRhoVals( s, p, xVals, fVals, dVals )
	yVals = abs(xVals-s).^p;
	if (~isrealarray(yVals,size(xVals)))
		rhoVals = -1.0;
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
		errFlag = true;
		return;
		% For example, 47^209 = "Inf".
	end
	vecC = matB\(matA'*matD*(fVals'));
	rhoVals = vecC(1) + (vecC(2)*yVals) - fVals;
	errFlag = false;
return;
end
