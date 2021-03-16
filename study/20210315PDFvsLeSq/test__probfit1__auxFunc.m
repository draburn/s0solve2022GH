clear;
setprngstates(0);
%funchP0 = @(c0,c1,cr)( ones(size(c0)) );
funchP0 = @(c0,c1,cr)( exp(-(cr.^2)) );
funchPr = @(r)( exp(-(r.^2)) );
xVals = [1,2,3,4];
%yVals = [4,3,2,1];
noiseLevel = 1E-1
yVals = 5-xVals + noiseLevel*randn(size(xVals));
wVals = ones(size(xVals));
%
c0Lo = 4;
c0Hi = 6;
c1Lo = -1.5;
c1Hi = -0.5;
crLo = noiseLevel/10.0;
crHi = noiseLevel*10.0;
%
if (0)
	c0Vals = linspace(c0Lo,c0Hi,11);
	c1Vals = linspace(c1Lo,c1Hi,11);
	cRVals = linspace(crLo,crHi,11);
	[ tensorC0, tensorC1, tensorCR ] = meshgrid(c0Vals,c1Vals,cRVals);
	vecC0 = reshape(tensorC0,1,[]);
	vecC1 = reshape(tensorC1,1,[]);
	vecCR = reshape(tensorCR,1,[]);
	tic();
	if (0)
	vecP = probfit1__auxFunc( ...
	  funchP0, funchPr, xVals, yVals, wVals, ...
	  vecC0, vecC1, vecCR );
	else
		assert(size(xVals,2)==4)
		funch_p = @(c0,c1,cr)(exp( -(cr.^2) ...
		  -(( yVals(1) - c0 - (c1*xVals(1)) )./cr).^2 ...
		  -(( yVals(2) - c0 - (c1*xVals(2)) )./cr).^2 ...
		  -(( yVals(3) - c0 - (c1*xVals(3)) )./cr).^2 ...
		  -(( yVals(4) - c0 - (c1*xVals(4)) )./cr).^2 ...
		  )./(cr.^4));
		vecP = funch_p(vecC0,vecC1,vecCR);
	end
	toc();
	%plot( vecC0, vecP, 'o' );
	%plot( vecC1, vecP, 'o' );
	plot( vecCR, vecP, 'o' );
	grid on;
	return;
end
%
tic();
if (0)
funch_p = @(c0,c1,cr)(probfit1__auxFunc( ...
  funchP0, funchPr, xVals, yVals, wVals, ...
  c0, c1, cr ));
else
	assert(size(xVals,2)==4)
	funch_p = @(c0,c1,cr)(exp( -(cr.^2) ...
	  -(( yVals(1) - c0 - (c1*xVals(1)) )./cr).^2 ...
	  -(( yVals(2) - c0 - (c1*xVals(2)) )./cr).^2 ...
	  -(( yVals(3) - c0 - (c1*xVals(3)) )./cr).^2 ...
	  -(( yVals(4) - c0 - (c1*xVals(4)) )./cr).^2 ...
	  )./(cr.^4));
end
%
funch_temp = @(c0,c1,cr)( funch_p(c0,c1,cr) );
mu_temp = triplequad( funch_temp, c0Lo, c0Hi, c1Lo, c1Hi, crLo, crHi );
mu1 = mu_temp
toc();
%
funch_temp = @(c0,c1,cr)( c0 .* funch_p(c0,c1,cr) );
mu_temp = triplequad( funch_temp, c0Lo, c0Hi, c1Lo, c1Hi, crLo, crHi );
muc0 = mu_temp
toc();
%
funch_temp = @(c0,c1,cr)( c1 .* funch_p(c0,c1,cr) );
mu_temp = triplequad( funch_temp, c0Lo, c0Hi, c1Lo, c1Hi, crLo, crHi );
muc1 = mu_temp
toc();
%
%funch_temp = @(c0,c1,cr)( (c0.^2) .* funch_p(c0,c1,cr) );
%mu_temp = triplequad( funch_temp, c0Lo, c0Hi, c1Lo, c1Hi, crLo, crHi );
%muc0sq = mu_temp
%toc();
%
funch_temp = @(c0,c1,cr)( (c1.^2) .* funch_p(c0,c1,cr) );
mu_temp = triplequad( funch_temp, c0Lo, c0Hi, c1Lo, c1Hi, crLo, crHi );
muc1sq = mu_temp
toc();
%
funch_temp = @(c0,c1,cr)( c0.*c1 .* funch_p(c0,c1,cr) );
mu_temp = triplequad( funch_temp, c0Lo, c0Hi, c1Lo, c1Hi, crLo, crHi );
muc0c1 = mu_temp
toc();
%
c0_avg = muc0 / mu1
c1_avg = muc1 / mu1
%c0sq_avg = muc0sq / mu1
c1sq_avg = muc1sq / mu1
c0c1_avg = muc0c1 / mu1
xPlus = -c0c1_avg / c1sq_avg


% Output 20210315.2200
%noiseLevel =  0.10000
%mu1 =  3.6050
%Elapsed time is 160.722 seconds.
%muc0 =  17.647
%Elapsed time is 356.882 seconds.
%muc1 = -3.4661
%Elapsed time is 552.883 seconds.
%muc1sq =  3.3625
%Elapsed time is 752.045 seconds.
%muc0c1 = -17.039
%Elapsed time is 950.637 seconds.
%c0_avg =  4.8950
%c1_avg = -0.96146
%c1sq_avg =  0.93273
%c0c1_avg = -4.7264
%xPlus =  5.0673
%octave:145> format long
%octave:146> xVals
%xVals =
%
%   1   2   3   4
%
%octave:147> yVals
%yVals =
%
%   3.87751634726318   3.07638376124291   1.95809767766837   1.05046192120313
%
%octave:148> wVals
%wVals =
%
%   1   1   1   1
%ctave:149> c0_avg 
%c0_avg =  4.89499128897781
%octave:150> c1_avg 
%c1_avg = -0.961463967659236
%octave:151> c1sq_avg 
%c1sq_avg =  0.932729574793570
%octave:152> c0c1_avg 
%c0c1_avg = -4.72640674432975
%octave:153> xPlus 
%xPlus =  5.06728517252794
