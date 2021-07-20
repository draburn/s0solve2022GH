function [ deltaX, deltaP ] = extFit__getDelta( ...
  bigX, bigP, xVals, fVals, wVals = [], prm = [] )
	if (isempty(wVals))
		wVals = ones(size(xVals));
	end
	%
	epsX = mygetfield(  prm,  "epsX",  min([ ...
	  sqrt(eps)*(max(xVals)-min(xVals)), ...
	  sqrt(sqrt(eps))*min(abs(diff(xVals))) ])  );
	epsP = mygetfield( prm, "epsP", sqrt(eps) );
	%
	rVals0 = extFit__getRVals( bigX, bigP, xVals, fVals, wVals );
	rValsPX = extFit__getRVals( bigX+epsX, bigP, xVals, fVals, wVals );
	rValsMX = extFit__getRVals( bigX-epsX, bigP, xVals, fVals, wVals );
	rValsPP = extFit__getRVals( bigX, bigP+epsP, xVals, fVals, wVals );
	rValsMP = extFit__getRVals( bigX, bigP-epsP, xVals, fVals, wVals );
	%
	dRdXVals = 0.5 * ( rValsPX - rValsMX ) / epsX;
	dRdPVals = 0.5 * ( rValsPP - rValsMP ) / epsP;
	%
	sumRXSq = sum( wVals .* dRdXVals .* dRdXVals );
	sumRPSq = sum( wVals .* dRdPVals .* dRdPVals );
	sumRXRP = sum( wVals .* dRdPVals .* dRdXVals );
	sumRXF  = sum( wVals .* dRdXVals .* rVals0 );
	sumRPF  = sum( wVals .* dRdPVals .* rVals0 );
	%
	matH = [ sumRXSq, sumRXRP; sumRXRP, sumRPSq ];
	matD = diag(abs(diag(matH))); % abs() may be superfluous.
	vecG = -[ sumRXF; sumRPF ];
	mu = mygetfield( prm, "mu", 0.0 );
	vecDelta = (matH+mu*matD)\vecG;
	%
	deltaX = vecDelta(1);
	deltaP = vecDelta(2);
return;
end


%!test
%!	thisFile = "test extFit__getRVals"
%!	setprngstates();
%!	bigA = randn();
%!	bigB = randn();
%!	bigX = randn()*exp(3*randn())
%!	bigP = 0.3 + abs(randn())
%!	numPts = 4 + round(abs(randn()*exp(3*randn())));
%!	funchF = @(x)( bigA + bigB * abs(x-bigX).^bigP );
%!	xVals = randn(1,numPts).*exp(randn()) + bigX;
%!	fVals = funchF(xVals);
%!	bigX_guess = bigX+0.1*randn()
%!	bigP_guess = bigP+0.1*randn()
%!	[ deltaX, deltaP ] = extFit__getDelta( ...
%!	  bigX_guess, bigP_guess, xVals, fVals )
%!	"..."
%!	resX_before = bigX_guess - bigX
%!	resP_before = bigP_guess - bigP
%!	resX = bigX_guess + deltaX - bigX
%!	resP = bigP_guess + deltaP - bigP
%!	msg( thisFile, __LINE__, "CHECK ABOVE VALUES!" );
