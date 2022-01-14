function [ rhoVals, errFlag ] = extFit__calcRhoVals( s, p, xVals, fVals, dVals )
	yVals = abs(xVals-s).^p;
	if (~isrealarray(yVals,size(xVals)))
		rhoVals = [];
		errFlag = true;
		return;
		% For example, 47^209 = "Inf".
	end
	matY = [ ones(size(yVals))', yVals' ];
	matD = diag(dVals);
	matA = matD*matY;
	r = rcond(matA'*matA);
	if ( ~isrealscalar(r) || eps >= r )
		rhoVals = [];
		errFlag = true;
		return;
		% For example, 47^209 = "Inf".
	end
	vecC = matA\(matD*(fVals'));
	rhoVals = vecC(1) + (vecC(2)*yVals) - fVals;
	errFlag = false;
return;
end
