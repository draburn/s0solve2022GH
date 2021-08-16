function [ errFlag, rhoVals, omega, fExt, f1 ] = findExt_sym__calcRes( xExt, p, xVals, fVals, dVals )
	yVals = abs(xVals-xExt).^p;
	if (~isrealarray(yVals,size(xVals)))
		errFlag = true;
		rhoVals = [];
		omega = [];
		fExt = [];
		f1 = [];
		return;
		% For example, 47^209 = "Inf".
	end
	matY = [ ones(size(yVals))', yVals' ];
	matD = diag(dVals);
	matA = matD*matY;
	matB = matA'*matA;
	r = rcond(matB); % Instead, check if all yVals are nearly the same...? Slower?
	if ( ~isrealscalar(r) || eps >= r )
		errFlag = true;
		rhoVals = [];
		omega = [];
		fExt = [];
		f1 = [];
		return;
	end
	vecC = matB\(matA'*matD*(fVals'));
	fExt = vecC(1);
	f1 = vecC(2);
	rhoVals = fExt + (f1*yVals) - fVals;
	omega = 0.5*sum(( dVals.*rhoVals ).^2);
	errFlag = false;
return;
end
