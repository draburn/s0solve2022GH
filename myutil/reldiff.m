%  Function...
%    rd = reldiff( x, y, denom0 )
%  Calculates rd = ||x-y|| / sqrt( ||x||^2 + ||y||^2 + denom0 ) for arbitary shaped x and y.
%  If both x and y are exactly zero, rd = 0.0.

function rd = reldiff( x, y, denom0=0.0 )
	if ( 2>nargin || 3<nargin || 1<nargout )
		print_usage();
	end
	numer = sumsq(reshape( x - y, 1, [] ));
	denom = sumsq(reshape(x,1,[])) + sumsq(reshape(y,1,[])) + denom0;
	if ( 0.0 == denom )
		rd = 0.0;
		return;
	end
	rd = sqrt( numer / denom );
return;
end

%!test
%!	setprngstates();
%!	assert( 0.0 == reldiff(0.0,0.0) );
%!	%
%!	for numTrials=1:100
%!		numDim = 1;
%!		sz = [ 1 + round(3.0*abs(randn(1,numDim))), 1 ];
%!		x = randn(sz);
%!		y = randn(sz);
%!		assert(fleq(  reldiff(x,y),  sqrt( sum((x-y).^2)/sum((x.^2)+(y.^2)) )  ));
%!	end
%!	%
%!	for numTrials=1:100
%!		numDim = 2;
%!		sz = 1 + round(3.0*abs(randn(1,numDim)));
%!		x = randn(sz);
%!		y = randn(sz);
%!		assert(fleq(  reldiff(x,y),  sqrt( sum(sum((x-y).^2))/sum(sum((x.^2)+(y.^2))) )  ));
%!	end
%!	%
%!	for numTrials=1:100
%!		numDim = 3;
%!		sz = 1 + round(3.0*abs(randn(1,numDim)));
%!		x = randn(sz);
%!		y = randn(sz);
%!		assert(fleq(  reldiff(x,y),  sqrt( sum(sum(sum((x-y).^2)))/sum(sum(sum((x.^2)+(y.^2)))) )  ));
%!	end
%!	%
%!	for numTrials=1:100
%!		numDim = 3;
%!		sz = 1 + round(3.0*abs(randn(1,numDim)));
%!		x = randn(sz);
%!		y = randn(sz);
%!		denom0 = abs(randn());
%!		assert(fleq(  reldiff(x,y,denom0),  sqrt( sum(sum(sum((x-y).^2)))/(denom0+sum(sum(sum((x.^2)+(y.^2))))) )  ));
%!	end
%!	%
%!	for numTrials=1:100
%!		numDim = 1 + round(3.0*abs(randn()));
%!		sz = 1 + round(3.0*abs(randn(1,numDim)));
%!		x = randn(sz);
%!		y = randn(sz);
%!		reldiff(x,y); % Just make sure it runs.
%!	end
%!	%
%!	for numTrials=1:100
%!		numDim = 1 + round(3.0*abs(randn()));
%!		sz = 1 + round(3.0*abs(randn(1,numDim)));
%!		x = zeros(sz);
%!		y = zeros(sz);
%!		assert( 0.0 == reldiff(x,y) );
%!	end
