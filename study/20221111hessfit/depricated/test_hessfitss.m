clear
mydefs;
setprngstates( 0 );
%
sizeX = 100
numPts = 20
%numPts = 200
numUnk = (sizeX*(sizeX+1))/2 + sizeX + 1
numGvn = (1+sizeX)*numPts
%
secret_f0 = 0.0;
secret_vecX0 = randn(sizeX,1);
foo = randn(sizeX,sizeX);
secret_matH = foo'*foo;
clear foo;
%
matX = randn(sizeX,numPts);
matX(:,1) = 0.0;
secret_matDX = matX - secret_vecX0;
rvecF = secret_f0 + sum( secret_matDX .* ( secret_matH * secret_matDX ), 1 )/2.0;
matG = secret_matH * secret_matDX;
%
matG += randn(size(matG));

[ f_minWas, pt0 ] = min(rvecF);
%
prm = [];
prm.pt0 = pt0;
msg( __FILE__, __LINE__, "Calling hessfitss()..." );
tic();
[ matV, fss, vecGss, matHss, hessfitssDat ] = hessfitss( sizeX, numPts, matX, rvecF, matG, pt0, prm );
toc();
eigsOfMatHss = eig(matHss)
if ( min(eigsOfMatHss) <= 0.0 )
	%matHss +=  (abs(min(eigsOfMatHss)) + sqrt(eps)*max(abs(eigsOfMatHss)) )*eye(size(matHss)); %???
	foo = abs(min(eigsOfMatHss)) + 0.1*max(abs(eigsOfMatHss));
	matHss += foo*eye(size(matHss));
endif
f_minWas = f_minWas
vecYNewton = -(matHss\vecGss);
vecXNext = matX(:,pt0) + matV * vecYNewton;
secret_vecDXNext = vecXNext - secret_vecX0;
%
if (1)
	vecYTest = -0.001*vecGss;
	f_test = fss + ( vecGss'*vecYTest ) + ( vecYTest'*matHss*vecYTest )/2.0
	vecXTest = matX(:,pt0) + matV*vecYTest;
	secret_vecDXTest = vecXTest - secret_vecX0;
	f_testWouldBe = secret_f0 + (secret_vecDXTest'*secret_matH*secret_vecDXTest)/2.0
endif
%
f_expect = fss + ( vecGss'*vecYNewton ) + ( vecYNewton'*matHss*vecYNewton )/2.0
f_wouldBe = secret_f0 + (secret_vecDXNext'*secret_matH*secret_vecDXNext)/2.0
%
msg( __FILE__, __LINE__, "End of test." ); return;
