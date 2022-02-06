function testFuncPrm = testfunc2021_genPrm( ...
  sizeX = 2, ...
  sizeF = 2, ...
  prngstates = 0, ...
  forceBadMin = true, ...
  forceSingleExt = true, ...
  forceGoodBehavior = true )
	%
	assert(isposintscalar(sizeX));
	assert(isposintscalar(sizeF));
	assert(isbool(forceBadMin));
	assert(issize(forceBadMin,[1,1]));
	assert(isbool(forceSingleExt));
	assert(issize(forceSingleExt,[1,1]));
	assert(isbool(forceGoodBehavior));
	assert(issize(forceGoodBehavior,[1,1]));
	if (~isempty(prngstates))
		setprngstates(prngstates);
	else
		setprngstates();
	end
	%
	%
	%
	vecXE = randn(sizeX,1);
	vecFE = randn(sizeF,1);
	matJPre = randn(sizeF,sizeX);
	ary3K = randn(sizeX,sizeX,sizeF);
	%
	testFuncPrm.vecXE = vecXE;
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	if (forceBadMin)
		testFuncPrm.vecFE = vecFE;
		testFuncPrm.matJ = matJPre - vecFE*(vecFE'*matJPre)/(vecFE'*vecFE);
		vecFEHat = vecFE/norm(vecFE);
		%
		if (forceSingleExt)
			for n=1:sizeX
			for m=1:sizeX
			if ( n==m )
				ary3K(n,m,:) = norm(reshape(ary3K(n,m,:),1,[]))*vecFEHat;
			else
				ary3K(n,m,:) = 0.0;
			end
			end
			end
		elseif (forceGoodBehavior)
			for n=1:sizeX
			for m=1:sizeX
				vecT1 = ary3K(n,m,:);
				vecT1 = reshape(vecT1,sizeF,1);
				vecT2 = vecT1 - vecFEHat*(vecFEHat'*vecT1);
				if (n==m)
					ary3K(n,m,:) = vecT2 + norm(reshape(ary3K(n,m,:),1,[]))*vecFEHat;
				else
					ary3K(n,m,:) = vecT2;
				end
			end
			end
		end
		testFuncPrm.ary3K = ary3K;
	else
		testFuncPrm.vecFE = zeros(sizeF,1);
		testFuncPrm.matJ = matJPre;
		testFuncPrm.ary3K = ary3K;
	end
	%
	for n=1:sizeF
		matTemp = testFuncPrm.ary3K(:,:,n);
		testFuncPrm.ary3K(:,:,n) = ( matTemp' + matTemp ) / 2.0;
	end
return;
end
