%  Function...
%    xOut = cent( xIn )
%  Overview...
%    Part of myutil module.
%    Returns the algebraic average values between successive elements of xIn.
%    Intended to complement "diff(xIn)".
%  Input values...
%    xIn: Input vector of numerical values.
%  Output values...
%    xOut: Output vector of numerical values of the algebraic averages
%      between successive values of xIn.
function xOut = cent( xIn )
	xOut = ( xIn(1:end-1) + xIn(2:end) ) / 2.0;
return;
end

%!test
%!	x0 = 0.0;
%!	x1 = 1.0;
%!	sizeX = 10;
%!	vecX = linspace( x0, x1, sizeX );
%!	vecCX = cent( vecX );
%!	for n=1:sizeX-1
%!		xExpected = x0 + ((x1-x0)*(n-0.5)/(sizeX-1.0));
%!		xTol = (abs(vecCX(n)) + abs(xExpected))*sqrt(eps);
%!		%disp([ vecCX(n), xExpected, xTol ]);
%!		assert( abs( vecCX(n) - xExpected ) <= xTol )
%!	end
%!	plot( ...
%!	  (1:sizeX), vecX, 'o-', ...
%!	  (1:sizeX-1)+0.5, vecCX, 'x' );
%!	disp( "Does this graph look reasonable?" );
%!	grid on;
%
%!test
%!	x0 = randn;
%!	x1 = randn;
%!	sizeX = 2+ceil(abs(randn));
%!	vecX = linspace( x0, x1, sizeX );
%!	vecDX = cent( vecX );
%!	for n=1:sizeX-1
%!		xExpected = x0 + ((x1-x0)*(n-0.5)/(sizeX-1.0));
%!		xTol = (abs(vecDX(n)) + abs(xExpected))*sqrt(eps);
%!		%disp([ vecDX(n), xExpected, xTol ]);
%!		assert( abs( vecDX(n) - xExpected ) <= xTol )
%!	end
