ax = [];
sizeX = 2;
sizeF = 2;
caseNum = 901;
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
	subCaseNum = 20
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
	sizeF = 2;
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
	elseif (1)
		%setprngstates(61937600); % A case with a sepx.
		%setprngstates(57200560); % sepx
		setprngstates(80098080); % Yeah, okay.
		temp_calpha = 10000.0;
		temp_cbeta = 0.01;
	else
		%setprngstates(3672896); % sepx
		setprngstates(21322288); % sepx in correct place.
		temp_calpha = 3.0;
		temp_cbeta = 0.3;
	end
	%
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	temp_matJ = randn(sizeF,sizeX);
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
	temp_vecF0 = temp_calpha * ( temp_matU * (temp_matU'*temp_vecLambda) ) ...
	  - temp_cbeta * ( temp_vecVHat * (temp_vecVHat'*temp_vecLambda) );
	%
	vecX0 = randn(sizeX,1);
	testFuncPrm.vecXE = vecX0;
	testFuncPrm.vecFE = temp_vecF0;
otherwise
	error( "Invalid value of switch." );
end
