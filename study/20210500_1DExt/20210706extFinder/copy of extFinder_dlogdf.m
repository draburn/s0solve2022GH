function [ bigX, bigP, datOut ] = extFinder_dlogdf( xVals, fVals, prm=[], datIn=[] );
	thisFile = "extFinder_dlogdf";
	%
	extFinder__init;
	thisFile = "extFinder_dlogdf";
	if (haveCand)
		thisFile = [ "RETURN from " thisFile ];
		return;
	end
	%
	if ( haveLeft && haveRight )
		n0 = nOfPtWiseMin-1;
		n1 = nOfPtWiseMin+1;
	elseif ( haveLeft )
		n0 = nOfPtWiseMin-1;
		n1 = nOfPtWiseMin;
	elseif ( haveRight )
		n0 = nOfPtWiseMin;
		n1 = nOfPtWiseMin+1;
	else
		error( "IMPOSSIBLE: have neither left nor right." );
	end
	
	%n0 = 2
	%n1 = numPts-1
	
	x0 = xVals(n0);
	x1 = xVals(n1);
	p0 = polyfit( xVals(n0-1:n0+1), fVals(n0-1:n0+1), 2 );
	p1 = polyfit( xVals(n1-1:n1+1), fVals(n1-1:n1+1), 2 );
	%
	%
	df0 = 2.0*p0(1)*x0 + p0(2);
	ddf0 = 2.0*p0(1);
	if ( 0.0==df0 )
		h0 = 0.0;
	else
		h0 = df0./abs(sqrt(eps)*abs(df0)+abs(ddf0));
	end
	df1 = 2.0*p1(1)*x1 + p1(2);
	ddf1= 2.0*p1(1);
	if ( 0.0==df1 )
		h1 = 0.0;
	else
		h1 = df1./abs(sqrt(eps)*abs(df1)+abs(ddf1));
	end
	assert( h1 > h0 );
	%
	%
	bigX = ( h1*x0 - h0*x1 ) / ( h1 - h0 )
	bigP = 1.0 + ( x1 - x0 ) / ( h1 - h0 )
	%
	%
thisFile = [ "RETURN from " thisFile ];	
return;
end
