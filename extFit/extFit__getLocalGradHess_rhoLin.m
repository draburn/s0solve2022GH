function datOut = extFit__getLocalGradHess_rhoLin( ...
  bigX, bigP, xVals, fVals, wVals = [], prm = [] )
	thisFile = "extFit__getLocalGradHess_rhoLin";
	if (isempty(wVals))
		wVals = ones(size(xVals));
	end
	%
	epsX = mygetfield(  prm,  "epsX",  min([ ...
	  sqrt(eps)*(max(xVals)-min(xVals)), ...
	  sqrt(sqrt(eps))*min(abs(diff(xVals))) ])  );
	epsP = mygetfield( prm, "epsP", sqrt(eps) );
	%
	rhoVals0  = extFit__getRhoVals( bigX, bigP, xVals, fVals, wVals );
	rhoValsPX = extFit__getRhoVals( bigX+epsX, bigP, xVals, fVals, wVals );
	rhoValsMX = extFit__getRhoVals( bigX-epsX, bigP, xVals, fVals, wVals );
	rhoValsPP = extFit__getRhoVals( bigX, bigP+epsP, xVals, fVals, wVals );
	rhoValsMP = extFit__getRhoVals( bigX, bigP-epsP, xVals, fVals, wVals );
	%
	rhoDXVals = ( rhoValsPX - rhoValsMX ) / ( 2.0 * epsX );
	rhoDPVals = ( rhoValsPP - rhoValsMP ) / ( 2.0 * epsP );
	%
	sumRXSq = sum( wVals .* rhoDXVals .* rhoDXVals );
	sumRPSq = sum( wVals .* rhoDPVals .* rhoDPVals );
	sumRXRP = sum( wVals .* rhoDPVals .* rhoDXVals );
	sumRXF  = sum( wVals .* rhoDXVals .* rhoVals0 );
	sumRPF  = sum( wVals .* rhoDPVals .* rhoVals0 );
	%
	vecG = -[ sumRXF; sumRPF ];
	matH = [ sumRXSq, sumRXRP; sumRXRP, sumRPSq ];
	%
	% Copy main stuff to datOut...
	datOut.vecG = vecG;
	datOut.matH = matH;
	%
	% Copy other stuff to datOut...
	datOut.bigX = bigX;
	datOut.bigP = bigP;
	datOut.xVals = xVals;
	datOut.fVals = fVals;
	datOut.wVals = wVals;
	datOut.prm = prm;
	datOut.epsX = epsX;
	datOut.epsP = epsP;
	datOut.rhoVals0 = rhoVals0;
	datOut.rhoValsPX = rhoValsPX;
	datOut.rhoValsMX = rhoValsMX;
	datOut.rhoValsPP = rhoValsPP;
	datOut.rhoValsMP = rhoValsMP;
	datOut.rhoDXVals = rhoDXVals;
	datOut.rhoDPVals = rhoDPVals;
	datOut.sumRXSq = sumRXSq;
	datOut.sumRPSq = sumRPSq;
	datOut.sumRXRP = sumRXRP;
	datOut.sumRXF  = sumRXF;
	datOut.sumRPF  = sumRPF;
return;
end

%!test
%!	thisFile = "test extFit__getLocalGradHess_rhoLin"
%!	bigA = 1.0;
%!	bigB = 1.0;
%!	bigX = 0.5;
%!	bigP = 2.0;
%!	funchF = @(x)( bigA + bigB * abs(x-bigX).^bigP );
%!	xVals = linspace( -2.0, 3.0, 6 );
%!	fVals = funchF(xVals);
%!	bigX_guess = bigX + 0.1;
%!	bigP_guess = bigP + 0.1;
%!	datOut = extFit__getLocalGradHess_rhoLin( bigX_guess, bigP_guess, xVals, fVals )
%!	1.0 - (datOut.matH(1,2)^2)/(datOut.matH(1,1)*datOut.matH(2,2))
%!	msg( thisFile, __LINE__, "CHECK ABOVE VALUES!" );
