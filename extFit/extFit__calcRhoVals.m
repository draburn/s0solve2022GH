function [ rhoVals, exceptionFlag ] = extFit__calcRhoVals( s, p, xVals, fVals, dVals )
	yVals = abs(xVals-s).^p;
	matY = [ ones(size(yVals))', yVals' ];
	matD = diag(dVals);
	matA = matD*matY;
	matB = matA'*matA;
	if ( eps^3>=rcond(matB) )
		exceptionFlag = true;
		return;
		% For example, 47^209 = "Inf".
	end
	vecC = matB\(matA'*matD*(fVals'));
	rhoVals = vecC(1) + (vecC(2)*yVals) - fVals;
	exceptionFlag = false;
return;
end
