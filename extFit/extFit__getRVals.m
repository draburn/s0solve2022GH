function rVals = extFit__getRVals( bigX, bigP, xVals, fVals, wVals = [] )
	if (isempty(wVals))
		wVals = ones(size(xVals));
	end
	vecC = extFit__getCoeff( bigX, bigP, xVals, fVals, wVals );
	dVals = abs( xVals - bigX ).^bigP;
	rVals = vecC(1) + (vecC(2) * dVals) - fVals;
return;
end


%!test
%!	thisFile = "test extFit__getRVals"
%!	setprngstates();
%!	bigA = randn();
%!	bigB = randn();
%!	bigX = randn()*exp(3*randn());
%!	bigP = 0.3 + abs(randn());
%!	numPts = 2 + round(abs(randn()*exp(3*randn())));
%!	funchF = @(x)( bigA + bigB * abs(x-bigX).^bigP );
%!	xVals = randn(1,numPts).*exp(randn()) + bigX;
%!	fVals = funchF(xVals);
%!	noiseLevel = 0.01*abs(randn())*sum(abs(fVals))/numPts;
%!	fVals += noiseLevel*randn(1,numPts)
%!	rVals = extFit__getRVals( bigX, bigP, xVals, fVals )
%!
%!	msg( thisFile, __LINE__, "PLEASE CHECK THESE VALUES!" );
