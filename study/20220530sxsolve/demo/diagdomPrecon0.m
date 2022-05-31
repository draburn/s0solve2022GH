% See Notes 2022-05-30-2100.
clear;
numFigs = 0;
%
sizeX = 100;
sizeV = 3;
setprngstates(0);
%
sizeF = sizeX;
%matJ = eye(sizeX,sizeX);
%matJ = diag((1:sizeX));
matJ = diag((1:sizeX)) + 1.0e-1*randn(sizeF,sizeX);
matV = randn(sizeX,sizeV);
matV = orth(matV);
assert( reldiff(matV'*matV,eye(sizeV,sizeV)) < sqrt(eps) );
matW = matJ*matV;
%
vScale = sqrt(sum(sum(matV.^2)));
wScale = sqrt(sum(sum(matW.^2)));
c1 = 1.0/(vScale^2);
c2 = 1.0/(wScale^2);
alpha = vScale/wScale;
eps1 = sqrt(eps);
eps2 = sqrt(eps);
%
matW2 = matW.^2;
matV2 = matV.^2;
matVW = matV.*matW;
%
if (0)
	numFigs++; figure(numFigs);
	aVals = linspace( -2.0, 2.0, 1001 );
	plot( aVals, 0*aVals, 'k-' );
	hold on;
	grid on;
	for m=1:sizeX
		w2 = sum(matW2(m,:));
		vw = sum(matVW(m,:));
		v2 = sum(matV2(m,:));
		%
		fVals = (c1*w2+eps1/(alpha^2))*aVals - (c1*vw+eps1/alpha) + (c2*vw+eps2*alpha)./(aVals.^2) - (c2*v2+eps2*alpha^2)./(aVals.^3);
		plot( aVals, fVals, 'o-' );
	endfor
	axis([ min(aVals), max(aVals), -alpha, alpha ]);
	hold off;
endif
%
if (1)
	numFigs++; figure(numFigs);
	aVals = linspace( -2.0, 2.0, 1001 );
	plot( aVals, 0*aVals, 'k-' );
	hold on;
	grid on;
	for m=1:sizeX
		w2 = sum(matW2(m,:));
		vw = sum(matVW(m,:));
		v2 = sum(matV2(m,:));
		%
		f0Vals = (c1*w2+eps1/(alpha^2))*aVals - (c1*vw+eps1/alpha) + (c2*vw+eps2*alpha)./(aVals.^2) - (c2*v2+eps2*alpha^2)./(aVals.^3);
		fVals = f0Vals.*(aVals.^3);
		plot( aVals, fVals, 'o-' );
	endfor
	axis([ min(aVals), max(aVals), min(fVals), -min(fVals) ]);
	hold off;
endif
%
for m=1:sizeX
	w2 = sum(matW2(m,:));
	vw = sum(matVW(m,:));
	v2 = sum(matV2(m,:));
	%
	foo = realroots([ c1*w2+eps1/(alpha^2), -(c1*vw+eps1/alpha), 0.0, c2*vw+eps2*alpha, -(c2*v2+eps2*alpha^2) ]);
	if ( vw > sqrt(eps)*sqrt(w2*v2) )
		foo = foo(foo>0.0);
	elseif ( vw < -sqrt(eps)*sqrt(w2*v2) )
		foo = foo(foo<0.0);
	else
		error( "Possible value is too close to zero." );
	endif
	assert( length(foo) == 1 );
	a(m) = foo;
endfor
matA = diag(a);
rcondComparison = [ rcond(matJ), rcond(matA*matJ) ]
