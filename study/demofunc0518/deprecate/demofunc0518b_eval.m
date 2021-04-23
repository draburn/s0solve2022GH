% 2D, support multi arg, numer and denom.

function matF = demofunc0518b_eval( matX, funcPrm )
	matFNumer = demofunc0518_eval( matX, funcPrm.numerPrm );
	matFDenom = demofunc0518_eval( matX, funcPrm.denomPrm );
	matF = matFNumer ./ matFDenom;
end
