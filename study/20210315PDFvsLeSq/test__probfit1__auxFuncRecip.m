"SLOWER!"

clear;
setprngstates(0);
xVals = [1,2,3,4];
%yVals = [4,3,2,1];
noiseLevel = 1E-3
yVals = 5-xVals + noiseLevel*randn(size(xVals));
%wVals = ones(size(xVals));
%
c0Lo = 4;
c0Hi = 6;
c1Lo = -1.5;
c1Hi = -0.5;
srLo = 1e-2/noiseLevel;
srHi = 10.0/noiseLevel;
%
if (0)
	c0Vals = 5%linspace(c0Lo,c0Hi,11);
	c1Vals = -1%linspace(c1Lo,c1Hi,11);
	sRVals = linspace(srLo,srHi,101);
	[ tensorC0, tensorC1, tensorSR ] = meshgrid(c0Vals,c1Vals,sRVals);
	vecC0 = reshape(tensorC0,1,[]);
	vecC1 = reshape(tensorC1,1,[]);
	vecSR = reshape(tensorSR,1,[]);
	tic();
	if (0)
	vecP = probfit1__auxFunc( ...
	  funchP0, funchPr, xVals, yVals, wVals, ...
	  vecC0, vecC1, vecSR );
	else
		assert(size(xVals,2)==4)
		funch_p = @(c0,c1,sr)(exp( -(1.0./(sr.^2)) ...
		  -(( yVals(1) - c0 - (c1*xVals(1)) ).*sr).^2 ...
		  -(( yVals(2) - c0 - (c1*xVals(2)) ).*sr).^2 ...
		  -(( yVals(3) - c0 - (c1*xVals(3)) ).*sr).^2 ...
		  -(( yVals(4) - c0 - (c1*xVals(4)) ).*sr).^2 ...
		  ).*(sr.^4));
		vecP = funch_p(vecC0,vecC1,vecSR);
	end
	toc();
	%plot( vecC0, vecP, 'o' );
	%plot( vecC1, vecP, 'o' );
	plot( vecSR, vecP, 'o' );
	grid on;
	return;
end
%
tic();
if (0)
funch_p = @(c0,c1,cr)(probfit1__auxFunc( ...
  funchP0, funchPr, xVals, yVals, wVals, ...
  c0, c1, cr ));
elseif (0)
	assert(size(xVals,2)==4)
	funch_p = @(c0,c1,cr)(exp( -(cr.^2) ...
	  -(( yVals(1) - c0 - (c1*xVals(1)) )./cr).^2 ...
	  -(( yVals(2) - c0 - (c1*xVals(2)) )./cr).^2 ...
	  -(( yVals(3) - c0 - (c1*xVals(3)) )./cr).^2 ...
	  -(( yVals(4) - c0 - (c1*xVals(4)) )./cr).^2 ...
	  )./(cr.^4));
else
	funch_p = @(c0,c1,sr)(exp( -(1.0./(sr.^2)) ...
	  -(( yVals(1) - c0 - (c1*xVals(1)) ).*sr).^2 ...
	  -(( yVals(2) - c0 - (c1*xVals(2)) ).*sr).^2 ...
	  -(( yVals(3) - c0 - (c1*xVals(3)) ).*sr).^2 ...
	  -(( yVals(4) - c0 - (c1*xVals(4)) ).*sr).^2 ...
	  ).*(sr.^4));
end
%
funch_temp = @(c0,c1,sr)( funch_p(c0,c1,sr) );
mu_temp = triplequad( funch_temp, c0Lo, c0Hi, c1Lo, c1Hi, srLo, srHi );
mu1 = mu_temp
toc();
%
funch_temp = @(c0,c1,sr)( c0 .* funch_p(c0,c1,sr) );
mu_temp = triplequad( funch_temp, c0Lo, c0Hi, c1Lo, c1Hi, srLo, srHi );
muc0 = mu_temp
toc();
%
funch_temp = @(c0,c1,sr)( c1 .* funch_p(c0,c1,sr) );
mu_temp = triplequad( funch_temp, c0Lo, c0Hi, c1Lo, c1Hi, srLo, srHi );
muc1 = mu_temp
toc();
%
%funch_temp = @(c0,c1,sr)( (c0.^2) .* funch_p(c0,c1,sr) );
%mu_temp = triplequad( funch_temp, c0Lo, c0Hi, c1Lo, c1Hi, srLo, srHi );
%muc0sq = mu_temp
%toc();
%
funch_temp = @(c0,c1,sr)( (c1.^2) .* funch_p(c0,c1,sr) );
mu_temp = triplequad( funch_temp, c0Lo, c0Hi, c1Lo, c1Hi, srLo, srHi );
muc1sq = mu_temp
toc();
%
funch_temp = @(c0,c1,sr)( c0.*c1 .* funch_p(c0,c1,sr) );
mu_temp = triplequad( funch_temp, c0Lo, c0Hi, c1Lo, c1Hi, srLo, srHi );
muc0c1 = mu_temp
toc();
%
c0_avg = muc0 / mu1
c1_avg = muc1 / mu1
%c0sq_avg = muc0sq / mu1
c1sq_avg = muc1sq / mu1
c0c1_avg = muc0c1 / mu1
xPlus = -c0c1_avg / c1sq_avg
