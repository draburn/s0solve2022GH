ax = [];
sizeX = 2;
sizeF = 2;
caseNum = 10010;
msg( __FILE__, __LINE__, sprintf( "caseNum = %d.", caseNum ) );
switch (caseNum)
case -1
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,[]);
	vecX0 = randn(2,1);
case 0
	% Easy linear case.
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 0.0; 0.0 ];
	testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 1.0 ];
	testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 0.0, 0.0; 0.0, 0.0 ];
	vecX0 = [ 1.0; 1.0 ];
case 1
	% Linear case.
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 0.0; 0.0 ];
	testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 3.0 ];
	testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 0.0, 0.0; 0.0, 0.0 ];
	vecX0 = [ 0.01; 3.0 ];
case 2
	% One-Component Quadratic.
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 0.1; 0.0 ];
	testFuncPrm.matJ = [ 0.0, 0.0; 0.0, 1.0 ];
	%%%testFuncPrm.ary3K(:,:,1) = [ 0.01, 0.0; 0.0, 0.0 ];
	%%%testFuncPrm.ary3K(:,:,2) = 10*[ 0.05, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,1) = [ 0.1, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 1.0, 0.0; 0.0, 0.0 ];
	%vecX0 = [ 5.0; 0.0 ];
	vecX0 = [ 2.2; -1.3 ];
	%vecX0 = [ 0.02; -0.0002 ];
	subCaseNum = 0
	switch (subCaseNum)
	case 0
		% Do nothing.
	case 10
		vecX0 = [ -0.5; 0.2 ]
	case 20
		vecX0 = [ -0.08; 0.02 ]
		ax = [ -0.1, 0.1, -0.01, 0.03 ]
	otherwise
		error( "Ivalid subCaseNum." );
	end
case 3
	% Narrow valley.
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 0.0; 0.0 ];
	testFuncPrm.matJ = [ 0.0, 0.0; 0.0, 1.0 ];
	%%%testFuncPrm.ary3K(:,:,1) = [ 0.01, 0.0; 0.0, 0.0 ];
	%%%testFuncPrm.ary3K(:,:,2) = 10*[ 0.05, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,1) = [ 0.1, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 0.5, 0.0; 0.0, 0.0 ];
	%vecX0 = [ 5.0; 0.0 ];
	vecX0 = [ 2.2; -1.3 ];
	%ax = [ -1.0, 5.5, -2, 0.5 ];
case 4
	% Linear narrow valley case.
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 0.0; 0.0 ];
	%testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 100.0 ];
	testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 10.0 ];
	testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 0.0, 0.0; 0.0, 0.0 ];
	vecX0 = [ 1.0; 1.0 ];
case 5
	% Linear narrow valley case.
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 0.0; 0.0 ];
	%testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 100.0 ];
	testFuncPrm.matJ = [ 1.1, 0.0; 0.0, 13.1 ];
	testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 0.0, 0.0; 0.0, 0.0 ];
	vecX0 = [ 1.0; 1.0 ];
case 20
	% Simple linear case.
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 0.0; 0.0 ];
	testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 3.0 ];
	testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 0.0, 0.0; 0.0, 0.0 ];
	vecX0 = [ 1.0; 1.0];
case 21
	% Simple linear case.
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 1.0; 1.0 ];
	testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 0.0, 0.0; 0.0, 1.0 ];
	vecX0 = [ 1.0; 1.0];
case 22
	% Simple linear case.
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 1.0; 1.0 ];
	testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 0.0, 0.0; 0.0, 0.0 ];
	vecX0 = [ 1.0; 1.0];
case 101
	% Dbl sepx.
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,25094192,true,false,false);
	vecX0 = [ 4.0; -0.11742 ]; % Some go to UR, one goes to LL.
case 102
	% Curves go down and split.
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,71832832,true,false,false);
	vecX0 = [ 3.0; 3.0 ];
	%ax = [ 2.413, 2.415, -2.74, -2.73 ];
	%ax = [ 2.411, 2.416, -2.75, -2.71 ];
	subCaseNum = 20
	switch (subCaseNum)
	case 0
		% Do nothing.
	case 10
		% OCQ is worse than cnst H.
		vecX0 = [ 2.2; -1.4 ]
	case 20
		% OCQ is better than cnst H.
		vecX0 = [ 2.4; -1.5 ]
	case 30
		vecX0 = [ 2.4; -2.6 ]
	otherwise
		error( "Ivalid subCaseNum." );
	end
case 103
	% Sharp narrow valley.
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,61409680,true,false,false);
	vecX0 = [ 3.0; 3.0 ];
case 104
	% Crescent
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,80724832,true,false,false);
	vecX0 = [ 3.0; 3.0 ];
case 500
	% OCQ jump-O "happens" to be good.
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,89071616);
	vecX0 = [ 1.0; 1.6 ];
case 510
	% Cnst J is better than OCQ.
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,98839248);
	vecX0 = randn(2,1);
case 520
	% Dramatically bad swerve.
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,46328864);
	vecX0 = randn(2,1);
case 530
	% Cnst J great; OCQ okay.
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,44135872);
	vecX0 = randn(2,1);
case 901
	% Try to make OCQ "1 -> 3 -> 1".
	sizeX = 2;
	sizeF = sizeX;
	assert( sizeX == sizeF );
	if (0)
		temp_calpha = 100.0; %Maybe?
		temp_cbeta = 0.01; %Maybe?
		% These results were using testFuncPrm.matJ = randn(sizeF,sizeX) indp of "temp_matJ"!
		%setprngstates(43661504); % Suggestive.
		%setprngstates(56075664); % From x-pt ish.
		%setprngstates(2017104); % All along same line.
		%
		% Full F still doesn't equal OCQ F, but, meh.
		%setprngstates(60599952);
	elseif (0)
		%setprngstates(61937600); % A case with a sepx.
		%setprngstates(57200560); % sepx
		setprngstates(80098080); % Yeah, okay.
		temp_calpha = 10000.0;
		temp_cbeta = 0.01;
	elseif (0)
		%setprngstates(3672896); % sepx
		setprngstates(21322288); % sepx in correct place.
		temp_calpha = 3.0;
		temp_cbeta = 0.3;
	elseif (1)
		setprngstates();
		temp_calpha = 10.0;
		temp_cbeta = 0.1;
	end
	%
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	temp_matJ = randn(sizeF,sizeX)
	testFuncPrm.matJ = temp_matJ;
	%
	temp_matH = temp_matJ' * temp_matJ
	[ temp_matEigVec, temp_matEigVal ] = eig(temp_matH)
	[ temp_vecEigValOfAbsMin, temp_nOfAbsMin ] = min(abs(diag(temp_matEigVal)))
	temp_vecPhiHat = temp_matEigVec(:,temp_nOfAbsMin)
	temp_vecEta = randn(sizeF,1);
	for n=1:sizeF
		testFuncPrm.ary3K(n,:,:) = temp_vecEta(n) * ( temp_vecPhiHat * (temp_vecPhiHat') );
	end
	%%%
	%%%temp_vecFoo1 = temp_matJ'*temp_matJ*temp_vecPhiHat
	%%%temp_vecFoo2 = temp_vecPhiHat*(temp_vecPhiHat'*temp_vecFoo1)
	%%%return;
	%%%
	%
	temp_matPsi = orth( eye(sizeX,sizeX) - (temp_vecPhiHat*(temp_vecPhiHat')), sqrt(eps) )
	assert( isrealarray(temp_matPsi,[sizeX,sizeX-1]) );
	temp_matW = temp_matJ * temp_matPsi
	temp_matU = orth( temp_matW )
	assert( isrealarray(temp_matU(sizeF,sizeX-1)) );
	temp_vecVHat = orth( eye(sizeX,sizeX) - (temp_matU*(temp_matU')), sqrt(eps) )
	assert( isrealarray(temp_vecVHat(sizeF,1)) );
	temp_matU'*temp_vecVHat
	assert( sum(sum(sum(sumsq(temp_matU'*temp_vecVHat)))) < sqrt(eps) );
	%
	temp_vecLambda = testFuncPrm.matJ * temp_vecPhiHat;
	%%%echo__UUTL = temp_matU*(temp_matU'*temp_vecLambda)
	%%%echo__VVTL = temp_vecVHat*(temp_vecVHat'*temp_vecLambda)
	temp_vecF0 = temp_calpha * ( temp_matU * (temp_matU'*temp_vecLambda) ) ...
	  - temp_cbeta * ( temp_vecVHat * (temp_vecVHat'*temp_vecLambda) )
	%
	%%%foo_f0tvl = (temp_vecVHat'*temp_vecF0)*(temp_vecVHat'*temp_vecLambda)
	%%%foo_f0til = temp_vecF0'*temp_vecLambda
	%
	vecX0 = randn(sizeX,1);
	testFuncPrm.vecXE = vecX0;
	testFuncPrm.vecFE = temp_vecF0;
case 920
	% Try to make OCQ "1 -> 3 -> 1" without phi being an eigenvector of JT J.
	sizeX = 2;
	sizeF = sizeX;
	assert( sizeX == sizeF );
	if (1)
		%setprngstates(8635232); % The forkpitch.
		%setprngstates(70144272); % A pitchfork.
		%setprngstates(15273696); % The FOCQ?
		setprngstates(34661776); % Indp closed curve!!!
		temp_calpha = 100.0;
		temp_cbeta = 0.01;
	end
	%
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	temp_matJ = randn(sizeF,sizeX)
	testFuncPrm.matJ = temp_matJ;
	%
	%temp_matH = temp_matJ' * temp_matJ
	%[ temp_matEigVec, temp_matEigVal ] = eig(temp_matH)
	%[ temp_vecEigValOfAbsMin, temp_nOfAbsMin ] = min(abs(diag(temp_matEigVal)))
	%temp_vecPhiHat = temp_matEigVec(:,temp_nOfAbsMin)
	temp_vecPhiHat = randn(sizeX,1);
	temp_vecPhiHat /= norm(temp_vecPhiHat);
	temp_vecEta = randn(sizeF,1);
	for n=1:sizeF
		testFuncPrm.ary3K(n,:,:) = randn(sizeX,sizeX);
	end
	applyExplicitSym = true
	if (applyExplicitSym)
	for k=1:sizeF
	for m=1:sizeX
	for n=1:m-1
		testFuncPrm.ary3K(k,m,n) += testFuncPrm.ary3K(k,n,m);
		testFuncPrm.ary3K(k,m,n) /= 2.0;
		testFuncPrm.ary3K(k,n,m) = testFuncPrm.ary3K(k,m,n);
	end
	end
	end
	end
	%%%
	temp_vecFoo1 = temp_matJ'*temp_matJ*temp_vecPhiHat
	temp_vecFoo2 = temp_vecPhiHat*(temp_vecPhiHat'*temp_vecFoo1)
	%%%return;
	%%%
	%
	temp_matPsi = orth( eye(sizeX,sizeX) - (temp_vecPhiHat*(temp_vecPhiHat')), sqrt(eps) )
	assert( isrealarray(temp_matPsi,[sizeX,sizeX-1]) );
	temp_matW = temp_matJ * temp_matPsi
	temp_matU = orth( temp_matW )
	assert( isrealarray(temp_matU(sizeF,sizeX-1)) );
	temp_vecVHat = orth( eye(sizeX,sizeX) - (temp_matU*(temp_matU')), sqrt(eps) )
	assert( isrealarray(temp_vecVHat(sizeF,1)) );
	temp_matU'*temp_vecVHat
	assert( sum(sum(sum(sumsq(temp_matU'*temp_vecVHat)))) < sqrt(eps) );
	%
	temp_vecLambda = testFuncPrm.matJ * temp_vecPhiHat;
	%%%echo__UUTL = temp_matU*(temp_matU'*temp_vecLambda)
	%%%echo__VVTL = temp_vecVHat*(temp_vecVHat'*temp_vecLambda)
	temp_vecF0 = temp_calpha * ( temp_matU * (temp_matU'*temp_vecLambda) ) ...
	  - temp_cbeta * ( temp_vecVHat * (temp_vecVHat'*temp_vecLambda) )
	%
	%%%foo_f0tvl = (temp_vecVHat'*temp_vecF0)*(temp_vecVHat'*temp_vecLambda)
	%%%foo_f0til = temp_vecF0'*temp_vecLambda
	%
	vecX0 = randn(sizeX,1);
	testFuncPrm.vecXE = vecX0;
	testFuncPrm.vecFE = temp_vecF0;
case 930
	% DRaburn 2022.01.31:
	%  better prepared here.
	% Let's try vecLambda = 0...
	%setprngstates(0); % Non-jump split.
	%setprngstates(66607104); % Causes an error. RESOLVED
	%setprngstates(13741184); % ERROR. RESOLVED
	%setprngstates(89520416); % Produces a fork 2321.
	setprngstates(73034896);
	%
	sizeX = 2;
	sizeF = sizeX;
	assert( sizeX == sizeF );
	temp_matIX = eye(sizeX,sizeX);
	%
	temp_matJ = randn(sizeF,sizeX);
	temp_vecPhiHat = randn(sizeX,1);
	temp_vecPhiHat /= norm(temp_vecPhiHat)
	temp_matPhiPhiT = temp_vecPhiHat*(temp_vecPhiHat');
	temp_matJ = temp_matJ*( temp_matIX - temp_matPhiPhiT )
	%
	temp_matPsi = orth( temp_matIX - temp_matPhiPhiT, sqrt(eps) )
	sizePsi = size(temp_matPsi,2);
	assert( isrealarray(temp_matPsi,[sizeX,sizePsi]) );
	assert( reldiff(temp_matPsi'*temp_matPsi,eye(sizePsi,sizePsi)) < sqrt(eps) );
	temp_matW = temp_matJ * temp_matPsi;
	temp_matQ = orth( temp_matW, sqrt(eps) )
	assert( isrealarray(temp_matQ,[sizeF,sizePsi]) );
	assert( reldiff(temp_matQ'*temp_matQ,eye(sizePsi,sizePsi)) < sqrt(eps) );
	%temp_matR = temp_matRQ'*temp_matW
	%
	temp_vecVHat = orth( temp_matIX - temp_matQ*(temp_matQ'), sqrt(eps) )
	assert( isrealarray(temp_vecVHat,[sizeX,1]) );
	assert( fleq(temp_vecVHat'*temp_vecVHat,1.0) );
	assert( sum(sumsq(temp_vecVHat'*temp_matQ)) < sqrt(eps) );
	if (1)
		temp_alpha = -1.0;
		temp_vecBeta = 5.0;
	elseif (1)
		temp_alpha = 0.1;
		temp_vecBeta = 10.0*(1.0+abs(randn(sizePsi,1)));
	else
		temp_alpha = 10.0;
		temp_vecBeta = 0.1*1.0./(1.0+abs(randn(sizePsi,1)));
	end
	%
	beWrong = false;
	if (beWrong)
		temp_vecF0 = temp_vecVHat*temp_alpha + temp_matPsi*temp_vecBeta;
		temp_vecEta = temp_vecVHat*temp_alpha - temp_matPsi*temp_vecBeta;
	else
		temp_vecF0 = temp_vecVHat*temp_alpha + temp_matQ*temp_vecBeta;
		temp_vecEta = temp_vecVHat*temp_alpha - temp_matQ*temp_vecBeta;
	end
	%
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	%
	vecX0 = randn(sizeX,1);
	testFuncPrm.vecXE = vecX0;
	testFuncPrm.vecFE = temp_vecF0;
	testFuncPrm.matJ = temp_matJ;
	for n=1:sizeF
		testFuncPrm.ary3K(n,:,:) = temp_vecEta(n) * ( temp_vecPhiHat * (temp_vecPhiHat') );
	end
	%
	temp_vecF0'*temp_vecVHat*temp_vecVHat'*temp_vecEta
	temp_vecF0'*temp_matQ*temp_matQ'*temp_vecEta
case 940
	% Simple ex of 1->3->1.
	% REALIZATION:
	%  Had thought ary3K was F*X*X,
	%  but it was actually X*X*F with a factor of 2.
	sizeX = 2;
	sizeF = 2;
	vecX0 = [ 0.0; 0.0 ];
	temp_matJ = [ 1.0, 0.0; 0.0, 0.0 ];
	%
	temp_vecPhiHat = [ 0.0; 1.0 ];
	temp_matPsi = [ 1.0; 0.0 ];
	temp_matQ = [ 1.0; 0.0 ];
	temp_vecVHat = [ 0.0; 1.0 ];
	%
	temp_vecF0 = [ -2.0; 1.0 ]
	temp_vecEta = [ 2.0; 1.0 ]
	%
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	testFuncPrm.vecXE = vecX0;
	testFuncPrm.vecFE = temp_vecF0;
	testFuncPrm.matJ = temp_matJ;
	for n=1:sizeF
		testFuncPrm.ary3K(:,:,n) = 2.0*temp_vecEta(n) * ( temp_vecPhiHat * (temp_vecPhiHat') );
	end
case 941
	% Alt ex of 1->3->1.
	sizeX = 2;
	sizeF = 2;
	vecX0 = [ 0.0; 0.0 ];
	temp_matJ = [ 1.0, 0.0; 0.0, 0.1 ];
	%
	temp_vecPhiHat = [ 0.0; 1.0 ];
	temp_matPsi = [ 1.0; 0.0 ];
	temp_matQ = [ 1.0; 0.0 ];
	temp_vecVHat = [ 0.0; 1.0 ];
	%
	temp_vecF0 = [ -2.0; 1.0 ]
	temp_vecEta = [ 2.0; 1.0 ]
	%
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	testFuncPrm.vecXE = vecX0;
	testFuncPrm.vecFE = temp_vecF0;
	testFuncPrm.matJ = temp_matJ;
	for n=1:sizeF
		testFuncPrm.ary3K(:,:,n) = 2.0*temp_vecEta(n) * ( temp_vecPhiHat * (temp_vecPhiHat') );
	end
	%
case 950
	% Alt ex of 1->3->1.
	sizeX = 2;
	sizeF = 2;
	vecX0 = [ 0.0; 0.0 ];
	%
	temp_matJ = [ 0.7338877403168940, -0.6195981167712608; 0.0143551588197744, 1.8483570954250721 ];
	temp_vecF0 = [ 4.36785955790190; 1.72861458692108 ];
	temp_vecPhiHat = [ -0.991770000512950; -0.128032285313289 ];
	temp_vecEta = [ -0.781166501131071; -0.246713234246678 ];
	%
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	testFuncPrm.vecXE = vecX0;
	testFuncPrm.vecFE = temp_vecF0;
	testFuncPrm.matJ = temp_matJ;
	for n=1:sizeF
		testFuncPrm.ary3K(:,:,n) = 2.0*temp_vecEta(n) * ( temp_vecPhiHat * (temp_vecPhiHat') );
	end
	%
case 960
	% Alt ex of 1->3->1.
	sizeX = 2;
	sizeF = 2;
	vecX0 = [ 0.0; 0.0 ];
	%
	temp_matJ = [ 3.1, 0.0; 0.0, -3.0 ];
	temp_vecF0 = [ -100.0; 1.0 ];
	temp_vecPhiHat = [ 0.0; 1.0 ];
	temp_vecEta = [ 0.1; 1.0 ];
	%
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	testFuncPrm.vecXE = vecX0;
	testFuncPrm.vecFE = temp_vecF0;
	testFuncPrm.matJ = temp_matJ;
	for n=1:sizeF
		testFuncPrm.ary3K(:,:,n) = 2.0*temp_vecEta(n) * ( temp_vecPhiHat * (temp_vecPhiHat') );
	end
case 10010
	% Copy of case 940, deomnstrating "'phi' 1->3->1";
	sizeX = 2;
	sizeF = 2;
	vecX0 = [ 0.0; 0.0 ];
	temp_matJ = [ 1.0, 0.0; 0.0, 0.0 ];
	%
	temp_vecPhiHat = [ 0.0; 1.0 ];
	temp_matPsi = [ 1.0; 0.0 ];
	temp_matQ = [ 1.0; 0.0 ];
	temp_vecVHat = [ 0.0; 1.0 ];
	%
	temp_vecF0 = [ -2.0; 1.0 ];
	temp_vecEta = [ 2.0; 1.0 ];
	%
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	testFuncPrm.vecXE = vecX0;
	testFuncPrm.vecFE = temp_vecF0;
	testFuncPrm.matJ = temp_matJ;
	for n=1:sizeF
		testFuncPrm.ary3K(:,:,n) = 2.0*temp_vecEta(n) * ( temp_vecPhiHat * (temp_vecPhiHat') );
	end
	%
case 10011
	% Similar to case 941, deomnstrating "disconnected closed loop".
	sizeX = 2;
	sizeF = 2;
	vecX0 = [ 0.0; 0.0 ];
	temp_matJ = [ 1.0, 0.0; 0.0, 0.02 ];
	%
	temp_vecPhiHat = [ 0.0; 1.0 ];
	temp_matPsi = [ 1.0; 0.0 ];
	temp_matQ = [ 1.0; 0.0 ];
	temp_vecVHat = [ 0.0; 1.0 ];
	%
	temp_vecF0 = [ -2.0; 1.0 ];
	temp_vecEta = [ 2.0; 1.0 ];
	%
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	testFuncPrm.vecXE = vecX0;
	testFuncPrm.vecFE = temp_vecF0;
	testFuncPrm.matJ = temp_matJ;
	for n=1:sizeF
		testFuncPrm.ary3K(:,:,n) = 2.0*temp_vecEta(n) * ( temp_vecPhiHat * (temp_vecPhiHat') );
	end
	%
case 10020
	% Copy of case 960, demonstrating "'backwards connected' 1->3->1"
	%  -- also called "1->3->1 reconnect" --
	%  though there is an additional "->3" afterwards;
	% Note: I suspect this "->3" is unavoidable in 2D when the vector "phiHat"
	% is taken to be the eigenvector corresponding to the smallest eigenvalue of J'*J.
	% Note: This suspicion is based on numerous trials with random values and manual attempts,
	% but I have not analyzed the algebra in detal.
	%
	sizeX = 2;
	sizeF = 2;
	vecX0 = [ 0.0; 0.0 ];
	%
	temp_matJ = [ 3.1, 0.0; 0.0, -3.0 ];
	temp_vecF0 = [ -100.0; 1.0 ];
	temp_vecPhiHat = [ 0.0; 1.0 ];
	temp_vecEta = [ 0.1; 1.0 ];
	%
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	testFuncPrm.vecXE = vecX0;
	testFuncPrm.vecFE = temp_vecF0;
	testFuncPrm.matJ = temp_matJ;
	for n=1:sizeF
		testFuncPrm.ary3K(:,:,n) = 2.0*temp_vecEta(n) * ( temp_vecPhiHat * (temp_vecPhiHat') );
	end
otherwise
	error( "Invalid value of switch." );
end

doKSymCheck = true;
for k=1:testFuncPrm.sizeF
for m=1:testFuncPrm.sizeX
for n=1:testFuncPrm.sizeX
	%%%assert( fleq( testFuncPrm.ary3K(k,m,n), testFuncPrm.ary3K(k,n,m) ) );
	assert( fleq( testFuncPrm.ary3K(m,n,k), testFuncPrm.ary3K(n,m,k) ) );
end
end
end
