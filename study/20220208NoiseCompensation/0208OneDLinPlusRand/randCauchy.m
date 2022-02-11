% Function...
%  x = randCauchy( sz, x1, x0 )
% Generates values from a Cuachy / Lorentz / "witch of Agnesi" distribution:
%  x = x0 + x1.*tan( pi*(rand(sz)-0.5) );
% Defaults are x0 = 0.0, x1 = 1.0, sz = 1.

function x = randCauchy( sz, x1, x0 )
	switch (nargin)
	case 0
		x = tan( pi*(rand()-0.5) );
	case 1
		x = tan( pi*(rand(sz)-0.5) );
	case 2
		x = x1.*tan( pi*(rand(sz)-0.5) );
	case 3
		x = x0 + x1.*tan( pi*(rand(sz)-0.5) );
	otherwise
		msg( __FILE__, __LINE__, "Bad nargin." );
		print_usage();
	end
return;
end

%!test
%!	numPts = 20001;
%!	x = randCauchy( [numPts,1] );
%!	% Let's sort and throw out some outliers...
%!	%xs = sort(x);
%!	%xs = xs(ceil(0.1*numPts):floor(0.9*numPts));
%!	hist( asinh(reshape(x,[],1)), 30 );
%!	grid on;
%!	xlabel( "ashin(x)" );
%!	ylabel( "hits" );
