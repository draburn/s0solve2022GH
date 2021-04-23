clear;
setprngstates(0);
xVals = [1,2,3,4]
noiseLevel = 1E-2
yVals = 5.0 - xVals + noiseLevel*randn(size(xVals))
probfit1__set_funcP;
%
c0Lo = 4.5
c0Hi = 5.5
c1Lo = -1.5
c1Hi = 0.5
srLo = 0.5
srHi = 2.0/noiseLevel
%
if (0)
	c0Vals = linspace(c0Lo,c0Hi,101);
	c1Vals = -1%linspace(c1Lo,c1Hi,101);
	sRVals = 1%linspace(srLo,srHi,101);
	[ tensorC0, tensorC1, tensorSR ] = meshgrid(c0Vals,c1Vals,sRVals);
	vecC0 = reshape(tensorC0,1,[]);
	vecC1 = reshape(tensorC1,1,[]);
	vecSR = reshape(tensorSR,1,[]);
	tic();
	vecP = funch_p(vecC0,vecC1,vecSR);
	toc();
	plot( vecC0, vecP, 'o' );
	%plot( vecC1, vecP, 'o' );
	%plot( vecSR, vecP, 'o' );
	grid on;
	return;
end
%
tic();
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
funch_temp = @(c0,c1,sr)( (c0.^2) .* funch_p(c0,c1,sr) );
mu_temp = triplequad( funch_temp, c0Lo, c0Hi, c1Lo, c1Hi, srLo, srHi );
muc0sq = mu_temp
toc();
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
c0sq_avg = muc0sq / mu1
c1sq_avg = muc1sq / mu1
c0c1_avg = muc0c1 / mu1
xPlus = -c0c1_avg / c1sq_avg

