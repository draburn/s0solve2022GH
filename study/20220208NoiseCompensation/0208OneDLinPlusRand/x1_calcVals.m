clear;
numFigs = 0;
tic();

msg( __FILE__, __LINE__, "NOPE. THIS IS WRONG. ABANDONING." );
%
% f = a + b*x + c*r.
caseNum = 100;
%
switch (caseNum)
case 0
	a_true = 0.0;
	b_true = 0.0;
	c_true = 0.0;
	xPts = [ 0.0, 1.0, 2.0 ];
	fPts = [ 0.0, 6.0, 6.0 ];
	numPts = 3;
case 1
	a_true = 0.0;
	b_true = 0.0;
	c_true = 0.0;
	xPts = [ 0.0, 1.0, 1.0 ];
	fPts = [ 0.0, 0.0, 2.0 ];
	numPts = 3;
case 100
	setprngstates(0);
	a_true = randCauchy();
	b_true = randCauchy();
	c_true = randCauchy()/1.0;
	numPts = 10;
	xPts = randn([1,numPts]);
	fPts = a_true+(b_true*xPts)+(c_true*randn([1,numPts]));
otherwise
	error( "Invalid caseNum." );
end
%
% Note: randn() has avg 0.0 and var 1.0; this is assumed in calculations below.
msg( __FILE__, __LINE__, sprintf( "Performing test with caseNum = %d...", caseNum ) );
msg( __FILE__, __LINE__, sprintf( "  a = %0.3e, b = %0.3e, c = %0.3e, N = %d.", a_true, b_true, c_true, numPts ) );
%
numFigs++; figure(numFigs);
plot( xPts, fPts, 'o', 'markersize', 10, 'linewidth', 2 );
%plot( ...
%  [min(xPts),max(xPts)], [0.0,0.0], 'k-', ...
%  [0.0,0.0], [min(fPts),max(fPts)], 'k-' );
grid on;
xlabel( "x" );
ylabel( "f" );
title( sprintf("f vs x %0.3e, %0.3e, %0.3e",a_true,b_true,c_true) );
%
% Use brute-force summation...
numerS = 0.0;
denomS = 0.0;
%for l=3:numPts
%for m=2:l-1
%for n=1:m-1
for l=1:numPts
for m=1:numPts
for n=1:numPts
	numerS += ( fPts(l)*(xPts(m)-xPts(n)) + fPts(m)*(xPts(n)-xPts(l)) + fPts(n)*(xPts(l)-xPts(m)) )^2;
	denomS += ( (xPts(m)-xPts(n)).^2 + (xPts(n)-xPts(l)).^2 + (xPts(l)-xPts(m)).^2 );
end
end
end
cSqAvg = numerS/denomS;
%cSqAvg /= numPts; % Why?
msg( __FILE__, __LINE__, sprintf( "c^2: expectation is %0.3e; true is %0.3e.", cSqAvg, c_true^2 ) );
%
numerS = 0.0;
denomS = 0.0;
normS = 0.0;
for m=2:numPts
for n=1:m-1
	numerS += (fPts(m)-fPts(n))^2;
	denomS += (xPts(m)-xPts(n))^2;
	normS += 1.0;
end
end
numerS = (numerS/normS) - (2.0*cSqAvg);
denomS = (denomS/normS);
bSqAvg = numerS/denomS;
msg( __FILE__, __LINE__, sprintf( "b^2: expectation is %0.3e; true is %0.3e.", bSqAvg, b_true^2 ) );
%
fxAvg = sum(fPts.*xPts)/numPts;
fAvg = sum(fPts)/numPts;
xAvg = sum(xPts)/numPts;
xSqAvg = sum(xPts.^2)/numPts;
xVarSq = xSqAvg - (xAvg^2);
%
aAvg = ( (fAvg*xSqAvg) - (fxAvg*xAvg) ) / xVarSq;
msg( __FILE__, __LINE__, sprintf( "a: expectation is %0.3e; true is %0.3e.", aAvg, a_true ) );
bAvg = ( fxAvg - (fAvg*xAvg) ) / xVarSq;
msg( __FILE__, __LINE__, sprintf( "b: expectation is %0.3e; true is %0.3e.", bAvg, b_true ) );
%
bVarSq = bSqAvg - (bAvg^2);
msg( __FILE__, __LINE__, sprintf( "bVarSq is %0.3e.", bVarSq ) );
if ( bVarSq < 0.0 )
	msg( __FILE__, __LINE__, "WARNING: bVarSq is neagative!" );
	warning( "bVarSq is negative!" );
end
%
%
toc();
msg( __FILE__, __LINE__, "AGAIN: THIS IS WRONG. ABANDONING." );
