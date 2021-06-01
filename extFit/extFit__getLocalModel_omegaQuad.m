function datOut = extFit__getLocalModel_omegaQuad( ...
  bigX, bigP, xVals, fVals, wVals = [], prm = [] )
	thisFile = "extFit__getLocalModel_omegaQuad";
	if (isempty(wVals))
		wVals = ones(size(xVals));
	end
	%
	epsX = mygetfield(  prm,  "epsX",  min([ ...
	  sqrt(eps)*(max(xVals)-min(xVals)), ...
	  sqrt(sqrt(eps))*min(abs(diff(xVals))) ])  );
	epsP = mygetfield( prm, "epsP", sqrt(eps) );
	%
	omega = zeros(3,3);
	for m=-1:1
	for n=-1:1
		omega(m+2,n+2) = 0.5 * sum( wVals.* ...
		  (extFit__getRhoVals( bigX + m*epsX, bigP + n*epsP, xVals, fVals, wVals )).^2 );
	end
	end
	%
	omegaX1 = ( omega(3,2) - omega(1,2) ) / (2.0*epsX);
	omegaP1 = ( omega(2,3) - omega(2,1) ) / (2.0*epsP);
	omegaX2 = ( omega(3,2) + omega(1,2) - 2.0*omega(2,2) ) / (epsX^2);
	omegaP2 = ( omega(2,3) + omega(2,1) - 2.0*omega(2,2) ) / (epsP^2);
	omegaXP = ( omega(3,3) + omega(1,1) - omega(3,1) - omega(1,3) ) ...
	  / (4.0*epsX*epsP);
	%
	matH = [ omegaX2, omegaXP; omegaXP, omegaP2 ];
	vecG = -[ omegaX1; omegaP1 ];
	prm_funchDelta = mygetfield( prm, "prm_funchDelta", [] );
	dat_funchDelta = extFit__getLocalModel__getFunchDelta( vecG, matH, prm_funchDelta );
	%
	% Copy main stuff to datOut...
	datOut.vecG = vecG;
	datOut.matH = matH;
	datOut.dat_funchDelta = dat_funchDelta;
	%
	% Copy input to datOut...
	datOut.bigX = bigX;
	datOut.bigP = bigP;
	datOut.xVals = xVals;
	datOut.fVals = fVals;
	datOut.wVals = wVals;
	datOut.prm = prm;
	%
	% Copy other stuff to datOut...
	datOut.epsX = epsX;
	datOut.epsP = epsP;
	%
	datOut.omega = omega;
return;
end

%!test
%!	thisFile = "test extFit__getLocalModel_omegaQuad"
%!	bigA = 1.0;
%!	bigB = 1.0;
%!	bigX = 0.5;
%!	bigP = 2.0;
%!	funchF = @(x)( bigA + bigB * abs(x-bigX).^bigP );
%!	xVals = linspace( -2.0, 3.0, 6 );
%!	fVals = funchF(xVals);
%!	bigX_guess = bigX + 0.1;
%!	bigP_guess = bigP + 0.1;
%!	datOut = extFit__getLocalModel_omegaQuad( bigX_guess, bigP_guess, xVals, fVals )
%!	1.0 - (datOut.matH(1,2)^2)/(datOut.matH(1,1)*datOut.matH(2,2))
%!	msg( thisFile, __LINE__, "CHECK ABOVE VALUES!" );
