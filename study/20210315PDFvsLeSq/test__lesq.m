clear;
setprngstates(0);
funchP0 = @(c0,c1,cr)( ones(size(c0)) );
funchPr = @(r)( exp(-(r.^2)) );
xVals = [1,2,3,4];
%yVals = [4,3,2,1];
noiseLevel = 1E-1
yVals = 5-xVals + noiseLevel*randn(size(xVals));
wVals = ones(size(xVals));
%
dat = lesq( [xVals',yVals',wVals'], 2 );
%
c0 = dat.vecC(1)
c1 = dat.vecC(2)
xPlus = -c0 / c1
