function vecP = probfit1__auxFunc( ...
  funchP0, ...
  funchPr, ...
  xVals, ...
  yVals, ...
  wVals, ...
  vecC0, ...
  vecC1, ...
  vecCR )
	%
	numVals = size(xVals,2);
	vecP = funchP0( vecC0, vecC1, vecCR );
	for n=1:numVals
		vecRho = (yVals(n) - vecC0 - (vecC1*xVals(n)))./vecCR;
		vecP .*= ( funchPr(vecRho) ./ vecCR ).^(wVals(n));
	end
	%
return;
end

%!test
%!	setprngstates(0);
%!	funchP0 = @(c0,c1,cr)( ones(size(c0)) );
%!	funchPr = @(r)( exp(-(r.^2)) );
%!	numVals = 3;
%!	xVals = [1,2,3];
%!	yVals = [3,2,1];
%!	wVals = ones(numVals,1);
%!	%
%!	vecC0 = linspace(0,10,1000);
%!	vecC1 = linspace(-5,5,1000);
%!	vecCR = linspace(eps,3,1000);
%!	vecP = probfit1__auxFunc( ...
%!	  funchP0, funchPr, xVals, yVals, wVals, ...
%!	  vecC0, vecC1, vecCR );
%!	%
%!	funch_p = @(c0,c1,cr)(probfit1__auxFunc( ...
%!	  funchP0, funchPr, xVals, yVals, wVals, ...
%!	  c0, c1, cr ));
%!	%
%!	tic();
%!	funch_temp = @(c0,c1,cr)( funch_p(c0,c1,cr) );
%!	mu_temp = triplequad( funch_temp, 1, 10, 1, 10, 1, 10 )
%!	toc();
%!	%
%!	funch_temp = @(c0,c1,cr)( c0 .* funch_p(c0,c1,cr) );
%!	mu_temp = triplequad( funch_temp, 1, 10, 1, 10, 1, 10 );
%!	toc();
%!	c0_avg = muc0 / mu1
