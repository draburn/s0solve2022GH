function rhoVals = extFit__getRhoVals( bigX, bigP, xVals, fVals, wVals = [] )
	thisFile = "extFit__getRhoVals";
	msg( thisFile, __LINE__, "DEPRECATED." );
	if (isempty(wVals))
		wVals = ones(size(xVals));
	end
	dVals = abs( xVals - bigX ).^bigP;
	sigma1 = sum(wVals);
	sigmaD = sum(wVals.*dVals);
	sigmaDD = sum(wVals.*dVals.*dVals);
	sigmaF = sum(wVals.*fVals);
	sigmaFD = sum(wVals.*dVals.*fVals);
	vecC = [ sigma1, sigmaD; sigmaD, sigmaDD ] \ [ sigmaF; sigmaFD ];
	rhoVals = vecC(1) + (vecC(2) * dVals) - fVals;
return;
end

%!test
%!	thisFile = "test extFit__getRhoVals"
%!	setprngstates();
%!	bigA = randn();
%!	bigB = randn();
%!	bigX = randn()*exp(3*randn());
%!	bigP = 0.3 + abs(randn());
%!	numPts = 5 + round(abs(randn()*exp(3*randn())));
%!	funchF = @(x)( bigA + bigB * abs(x-bigX).^bigP );
%!	xVals = randn(1,numPts).*exp(randn()) + bigX;
%!	fVals_before = funchF(xVals)
%!	noiseLevel = 0%0.01*abs(randn())*sum(abs(fVals_before))/numPts;
%!	noiseVals = noiseLevel*randn(1,numPts)
%!	fVals = fVals_before + noiseVals;
%!	rVals = extFit__getRhoVals( bigX, bigP, xVals, fVals )
%!	msg( thisFile, __LINE__, "ARE THESE VALUES FOR REASONABLE?" );
