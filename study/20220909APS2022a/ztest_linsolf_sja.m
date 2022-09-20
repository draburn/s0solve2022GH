clear;
mydefs;
setprngstates(59972704);
numFigs = 0;
sizeX = 10
sizeF = sizeX
%sizeK = ceil(sizeX^0.8)
sizeL = 2%ceil(sizeK/2.0)
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
%matV = randn(sizeX,sizeK);
%matV = utorthdrop( matV );
matV = [
   2.85286655819914e-01  -1.48139768928509e-01  -8.84754474687839e-02  -5.17527112723658e-03
   5.87862786241981e-01   3.53688180890744e-01   1.41213464792903e-01   1.02397993285685e-01
   5.24500180841158e-01  -5.91477970140600e-01  -2.64625838862542e-01  -2.17634867715639e-01
   4.53044920596477e-03   1.68410357206601e-03  -2.66009025896064e-03  -2.20687318568133e-03
   1.61112358310032e-01  -3.10304657641691e-01  -1.22381108010614e-01  -1.50731875209891e-01
   1.07569014516241e-02   9.60049436123856e-03   1.51373031132336e-02   3.51717994489680e-02
  -8.58075263680506e-05   4.24169034502771e-01  -5.90948293588992e-01  -5.83606134759250e-01
   1.79791640837740e-01   1.64237171342479e-01   6.44493504244832e-01  -6.46200258199894e-01
   4.76988457173242e-01   2.62825581966585e-01   4.74274903554469e-02   3.84553629438786e-01
  -1.09507859110022e-01  -3.61605585313260e-01   3.46625934961746e-01  -1.10038719473695e-01
];
assert( reldiff( matV'*matV, eye(4) ) < sqrt(eps) );
matW = matJSecret*matV;
%
prm = [];
time0 = time();
%[ matJA_basic, datOut ] = sja_scratch350( matV, matW, prm );
[ matJA_basic, datOut ] = sja_basic( matV, matW, prm );
time_basic= time()-time0
rdWV = reldiff( matJA_basic*matV, matW )
rd_basic = reldiff( matJSecret, matJA_basic )
%msg( __FILE__, __LINE__, sprintf( " sja_scratch100() produced %0.3e in %0.3f seconds.", rd_basic, time_basic ) );
%
funchF = @(x)( matJSecret*(x-vecXSecret) );
vecX0 = zeros(sizeX,1);
grootPrm.verbLev = VERBLEV__PROGRESS;
groot_jfnk_sja( funchF, vecX0, grootPrm );
%
groot_jfnk_baseline( funchF, vecX0, grootPrm );
