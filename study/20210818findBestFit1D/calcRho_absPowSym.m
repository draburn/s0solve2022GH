function [ errFlag, vecRho ] = calcRho_absPowSym( rhoArgs, vecZ )
	yVals = abs(rhoArgs.xVals-vecZ(1)).^vecZ(2);
	if (~isrealarray(yVals))
		errFlag = true;
		vecRho = [];
		return;
		% For example, 47^209 = "Inf".
	end
	matY = [ ones(size(yVals)), yVals ];
	matD = diag(rhoArgs.dVals);
	matA = matD*matY;
	matB = matA'*matA;
	r = rcond(matB); % Instead, check if all yVals are nearly the same...? Slower?
	if ( ~isrealscalar(r) || eps >= r )
		errFlag = true;
		vecRho = [];
		return;
	end
	vecC = matB\(matA'*matD*(rhoArgs.fVals));
	vecRho = (rhoArgs.dVals.^2) .*( vecC(1) + (vecC(2)*yVals) - rhoArgs.fVals );
	errFlag = false;
return;
end
%fExt = vecC(1);
%f1 = vecC(2);

%!test
%!	rhoArgs.xVals = [ 0.0, 1.0, 2.0, 3.0 ];
%!	rhoArgs.fVals = rhoArgs.xVals.^2;
%!	rhoArgs.dVals = ones(size(rhoArgs.xVals));
%!	%
%!	vecZ = [ 0.0, 0.0 ]
%!	[ errFlag, vecRho ] = calcRho_absPowSym( rhoArgs, vecZ );
%!	echo__errFlag = errFlag
%!	echo__vecRho = vecRho
%!	%
%!	vecZ = [ 1.0, 1.0 ]
%!	[ errFlag, vecRho ] = calcRho_absPowSym( rhoArgs, vecZ );
%!	echo__errFlag = errFlag
%!	echo__vecRho = vecRho
%!	%
%!	vecZ = [ 0.0, 2.0 ]
%!	[ errFlag, vecRho ] = calcRho_absPowSym( rhoArgs, vecZ );
%!	echo__errFlag = errFlag
%!	echo__vecRho = vecRho