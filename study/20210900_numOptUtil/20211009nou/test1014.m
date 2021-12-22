clear;
thisFile = "study1014"
tic();
setprngstates(53538896);
%%%setprngstates();
numFigs = 0;
%%%probSize = 100;
probSize = 5
for nTrials=1:1
matA = randn(probSize);
matB = randn(probSize);
matH = matA'*matA - 0.5*(matB'*matB);
%%%matH = matA'*matA - 0.001*(matB'*matB);
%matH = diag(abs(randn(probSize,1)));
%%%[ matPsi, matLambda ] = eig(matH)
vecG = randn(probSize,1);
%
vecDelta = zeros(probSize,1);
vecU = vecG + matH*vecDelta;
assert( norm(vecU) > 0.0 );
vecU /= norm(vecU);
matV = vecU;
clear vecU;
%
[ matR, cholFlag ] = chol( matV'*matH*matV );
%%%matHV = matH*matV;
%%%matVTH2V = matHV'*matHV;
%%%[ matR, cholFlag ] = chol( matVTH2V );
n = 0;
while (~cholFlag)
	n++;
	if ( n>= probSize )
		break;
	end
	%
	switch n
	case {1}
		h1 = matV'*matH*matV;
	case {2}
		h2 = matV'*matH*matV;
	case {3}
		h3 = matV'*matH*matV;
	case {4}
		h4 = matV'*matH*matV;
	otherwise
	end
	%
	%matV'*vecG % Just norm(vecG) * e1, I know.
	vecB = matV'*vecG;
	vecDelta = -matV * ( matR \ (matR'\(matV'*vecG)) );
	%
	% Use Krylov or not?
	%omega = vecG'*vecDelta + 0.5*vecDelta'*matH*vecDelta
	vecU = vecG + matH*vecDelta;
	%%%vecU = randn(probSize,1);
	normU = norm(vecU);
	assert( normU > 0.0 );
	vecU /= normU;
	for npass=1:2
		vecU -= matV * ( matV'*vecU );
		normU = norm(vecU);
		assert( normU > 0.0 );
		vecU /= normU;
	end
	matV = [ matV, vecU ];
	clear vecU;
	clear res;
	%
	[ matR, cholFlag ] = chol( matV'*matH*matV );
end
%[ matPsiVTHV, matLambdaVTHV ] = eig(matV'*matH*matV)
matV_original = matV;
if ( n>=probSize )
	%msg( thisFile, __LINE__, "Reached probSize." );
	continue;
end
%echo__n = n
%
% Swtich to other notation...
vecW = matV(:,end);
matV = matV(:,1:end-1);
if ( vecW'*matH*vecW < 0.0 )
	%msg( thisFile, __LINE__, "vecW'*matH*vecW < 0.0 is great, but not the case being examined here." );
	%msg( thisFile, __LINE__, "Alternate okay." );
	continue
end
%
%xi = vecW'*matH*vecW - vecW'*matH*matV*((matV'*matH*matV)\(matV'*matH*vecW))
%assert( xi >= 0.0 );
%
vecPsi = vecW - matV * ((matV'*matH*matV)\(matV'*matH*vecW));
assert( norm(vecPsi) > 0.0 );
vecPsi /= norm(vecPsi);
muCrit = -vecPsi'*matH*vecPsi;
%echo__muCrit = muCrit
assert( muCrit >= 0.0 );
%
matUAug = myorth( [vecPsi, matV] );
matU = matUAug(:,2:end);
matUTHU = matU'*matH*matU;
[ matPsiUTHU, matLambdaUTHU ] = eig(matUTHU);
%
muCritU = -min(min(matLambdaUTHU));
if ( muCritU < 0.0 )
	msg( thisFile, __LINE__, sprintf( "Good (%e, %e).", muCrit, muCritU ) );
elseif ( muCritU < eps )
	msg( thisFile, __LINE__, sprintf( "- Okay (%e, %e). -", muCrit, muCritU ) );
elseif ( muCritU <= muCrit )
	msg( thisFile, __LINE__, sprintf( "~~~ Meh (%e, %e). ~~~", muCrit, muCritU ) );
elseif ( muCritU > muCrit )
	msg( thisFile, __LINE__, sprintf( "***** BAD (%e, %e). *****", muCrit, muCritU ) );
	assert(0);
end
%
if (0)
[ matR, cholFlag ] = chol(matUTHU);
if (cholFlag)
	msg( thisFile, __LINE__, "matUTHU is NOT pos-def!!!!" );
	assert(0);
else
	msg( thisFile, __LINE__, "Okay." );
end
end
%
end %trials loop
echo__nTrials = nTrials
