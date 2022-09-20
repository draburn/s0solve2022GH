clear;
mydefs;
setprngstates(25279904);
numFigs = 0;
sizeX = 100
sizeF = sizeX
sizeK = ceil(sizeX^0.8)
sizeL = 5%ceil(sizeK/2.0)
%
vecXSecret = randn(sizeX,1);
matJSecret = eye(sizeF,sizeX);
for nf=1:sizeF
for m=1:sizeL-1
	nx = 1 + floor(sizeX*(1.0-100.0*eps)*rand);
	nzeListSecret(m) = nx;
	matJSecret(nf,nx) += randn();
endfor
endfor
%
matV = randn(sizeX,sizeK);
matV = utorthdrop( matV );
matW = matJSecret*matV;
%
if (0)
prm = [];
time0 = time();
[ matJA_basic, datOut ] = sja_scratch350( matV, matW, prm );
%[ matJA_basic, datOut ] = sja_basic( matV, matW, prm );
time_basic= time()-time0
rdWV = reldiff( matJA_basic*matV, matW )
rd_basic = reldiff( matJSecret, matJA_basic )
%msg( __FILE__, __LINE__, sprintf( " sja_scratch100() produced %0.3e in %0.3f seconds.", rd_basic, time_basic ) );
endif
%
funchF = @(x)( matJSecret*(x-vecXSecret) );
vecX0 = zeros(sizeX,1);
grootPrm.verbLev = VERBLEV__PROGRESS;
%
groot_jfnk_sja( funchF, vecX0, grootPrm );
groot_jfnk_baseline( funchF, vecX0, grootPrm );
