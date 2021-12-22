thisFile = "xinit";
%
ax = [];
sizeX = 2;
sizeF = 2;
caseNum = 103;
msg( thisFile, __LINE__, sprintf( "caseNum = %d.", caseNum ) );
switch (caseNum)
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
case 101
	% Dbl sepx.
	testFuncPrm = testFunc_genPrm(sizeX,sizeF,25094192,true,false,false);
	vecX0 = [ 4.0; -0.11742 ]; % Some go to UR, one goes to LL.
case 102
	% Curves go down and split.
	testFuncPrm = testFunc_genPrm(sizeX,sizeF,71832832,true,false,false);
	vecX0 = [ 3.0; 3.0 ];
	%ax = [ 2.413, 2.415, -2.74, -2.73 ];
	%ax = [ 2.411, 2.416, -2.75, -2.71 ];
case 103
	% Sharp narrow valley.
	testFuncPrm = testFunc_genPrm(sizeX,sizeF,61409680,true,false,false);
	vecX0 = [ 3.0; 3.0 ];
case 104
	% Crescent
	testFuncPrm = testFunc_genPrm(sizeX,sizeF,80724832,true,false,false);
	vecX0 = [ 3.0; 3.0 ];
otherwise
	error( "Invalid value of switch." );
end
%
thisFile = "RETURN FROM xinit";
return;
