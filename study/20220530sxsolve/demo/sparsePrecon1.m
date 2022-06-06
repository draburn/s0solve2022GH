% See Notes 2022-05-30-2100.
clear;
numFigs = 0;
setprngstates(0);
%
sizeF = 1;
sizeX = 1000;
numElemPerCol = 0;
numAddtlElemPerRow = 5;
c0 = 0.0;
csx = 0.0;
csf = 0.0;
%
matJ0 = zeros(sizeF,sizeX);
for n=1:sizeX
for t=1:numElemPerCol
	m = ceil( sqrt(eps) + (sizeF-2.0*sqrt(eps))*rand() );
	matJ0(m,n) = 1.0;
endfor
endfor
for m=1:sizeF
for t=1:numAddtlElemPerRow
	n = ceil( sqrt(eps) + (sizeX-2.0*sqrt(eps))*rand() );
	matJ0(m,n) = 1.0;
endfor
endfor
matJ0 += c0 * randn(sizeF,sizeX);
assert( isrealarray(matJ0,[sizeF,sizeX]) );

matSX = diag(exp(csx*randn(sizeX,1)));
matSF = diag(exp(csf*randn(sizeF,1)));
matJ = matSF*matJ0/matSX;



%setprngstates(3568384);
setprngstates(74504256);
sizeU0 = 50;

matU0 = randn(sizeX,sizeU0);
matV = utorthdrop(matU0);
if (0)
	matU3 = zeros(sizeX,3);
	matU3(1:3:end,1) = 1.0;
	matU3(2:3:end,2) = 1.0;
	matU3(3:3:end,3) = 1.0;
	matU5 = zeros(sizeX,5);
	matU5(1:5:end,1) = 1.0;
	matU5(2:5:end,2) = 1.0;
	matU5(3:5:end,3) = 1.0;
	matU5(4:5:end,4) = 1.0;
	matU5(5:5:end,5) = 1.0;
	matV = [ matV, matU3, matU5 ];
endif

%matV = matV(:,1:10);
%matV = matV(:,11:20);

sizeV = size(matV,2);
%assert( reldiff(matV'*matV,eye(sizeV,sizeV)) < sqrt(eps) );
matW = matJ*matV;


if (1)
	% GradObj doesn't help?!
	constraint_function = @(x)( matV'*x - matW(1,:)' );
	opts = optimset( "equc", {constraint_function}, "tolFun", 1e-16, "GradObj", "on" );
	objective_function = @(x)( getL1Norm(x) );
	pin = ( (matW(1,:)*(matV'))*pinv(matV*(matV')) )';
	[p, objf, cvg, outp] = nonlin_min (objective_function, pin, optimset ("equc", {constraint_function}, "tolFun", 1e-16 ));
	matJEst(1,:) = p;
else
objective_function = @(x)( norm(x,1) );
%%%objective_function = @(x)( sum(x.^2./(eps+abs(x))) );
pin = ( (matW(1,:)*(matV'))*pinv(matV*(matV')) )';
constraint_function = @(x)( matV'*x - matW(1,:)' );
[p, objf, cvg, outp] = nonlin_min (objective_function, pin, optimset ("equc", {constraint_function}, "tolFun", 1e-16 ));
for n=1:0
[p, objf, cvg, outp] = nonlin_min (objective_function, p, optimset ("equc", {constraint_function}, "tolFun", 1e-16 ));
endfor
matJEst(1,:) = p;
endif

objective_function( pin )
objective_function( p )
objective_function( matJ(1,:)' )
norm(constraint_function(pin))
norm(constraint_function(p))
norm(constraint_function(matJ(1,:)'))

matResCollective = matW - matJEst*matV;
if (1)
	numFigs++; figure(numFigs);
	plot( matJ(1,:), 'o-', 'linewidth', 3, p, 'x-' );
	grid on;
	numFigs++; figure(numFigs);
	plot( matResCollective, 'o-' );
	grid on;
endif

return;


 ## Example for default optimization (Levenberg/Marquardt with
 ## BFGS), one non-linear equality constraint. Constrained optimum is
 ## at p = [0; 1].
 objective_function = @ (p) p(1)^2 + p(2)^2;
 pin = [-2; 5];
 constraint_function = @ (p) p(1)^2 + 1 - p(2);
 [p, objf, cvg, outp] = nonlin_min (objective_function, pin, optimset ("equc", {constraint_function}))

return;
%
%
matJCollectiveEst = zeros(sizeF,sizeX);
matResCollective = zeros(sizeF,sizeV);
for m=1:sizeF
	funchOmega = @(x)( 1E8*norm(matV'*x-matW(m,:)') + norm(x,1) );
	%vecX0 = zeros(sizeX,1);
	%funchOmega(vecX0)
	vecX0 = matJ(m,:)'+0.1*randn(sizeX,1);
	%return
	fminunc_options = optimset( "TolFun", 1e-4 );
	[ vecXF, fminunc_fval, fminunc_info, fminunc_output, fminunc_grad, fminunc_hess ] = fminunc(funchOmega,vecX0,fminunc_options);
	%fminunc_fval
	fminunc_output
	fminunc_info
	funchOmega(vecX0)
	funchOmega(vecXF)
	funchOmega(matJ(m,:)')
	matJEst(m,:) = vecXF';
endfor
matResCollective = matW - matJEst*matV;
%matJCollectiveEst
%matJ

if (1)
	numFigs++; figure(numFigs);
	plot( matJ, 'o-', 'linewidth', 3, matJEst, 'x-' );
	grid on;
	numFigs++; figure(numFigs);
	plot( matResCollective, 'o-' );
	grid on;
endif
