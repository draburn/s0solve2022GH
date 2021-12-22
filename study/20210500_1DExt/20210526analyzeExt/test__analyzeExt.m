clear;
commondefs;
thisFile = "test__analyzeExt";
msg( thisFile, __LINE__, "THIS CODE IS UNWORTHY." );

setprngstates(57655840);
f0 = abs(randn());
f1 = abs(randn());
bigX = randn()*exp(abs(randn()));
funchF = @(x)( f0 + f1*abs(x-bigX).^6);
numPts = 100;
numTrials = 0;
while (1)
	numTrials++;
	foo1 = randn;
	foo2 = foo1*(1+abs(randn));
	foo3 = -foo1*(1+abs(randn));
	if ( 3==numPts )
		xVals = sort(bigX+[foo1,foo2,foo3]);
	else
		xVals = sort([ bigX+[foo1,foo2,foo3], randn(1,numPts-3) ]);
	end
	fVals = funchF(xVals);
	[ foo, indexOfFMin ] = min( fVals );
	if ( 1==indexOfFMin )
		continue;
	elseif ( numPts==indexOfFMin )
		continue;
	elseif ( fVals(indexOfFMin)<fVals(indexOfFMin-1) ...
	  && fVals(indexOfFMin)<fVals(indexOfFMin+1) )
		break;
	end
	assert( numTrials < 10 );
end

xVals = bigX+0.1+linspace(-2,2,100);
fVals = funchF(xVals);

numFigs = 0;
viz_numPts = 1000;
viz_xVals = linspace( ...
  min([min(xVals),bigX])-0.01, ...
  max([max(xVals),bigX])+0.01, ...
  viz_numPts );
viz_fVals = funchF(viz_xVals);

numFigs++; figure(numFigs);
plot( ...
  viz_xVals, viz_fVals, 'o-', 'markersize', 3, ...
  xVals, fVals, '+', 'linewidth', 5, 'markersize', 30, ...
  bigX, funchF(bigX), '*', 'linewidth', 4, 'markersize', 25  );
grid on;


[ xOfCand, meritOfCand, datOut ] = analyzeExt( xVals, fVals );
