%funchFOCQ_temp = @(dummyX)( vecF0 + matJ0*(dummyX-vecX0) + vecEta0_jtj*(vecPhi0_jtj'*(dummyX-vecX0))^2 );
funchFOCQ_temp = funchF;
funchOm_temp = @(dummyX)( sumsq(funchFOCQ_temp(dummyX))/2.0 );
vecX_temp = vecXG_jtj;
%vecX_temp = vecX0-2.0*vecG0/norm(vecG0)
%
vecF_temp = funchFOCQ_temp(vecX_temp);
%
matH_temp = zeros(sizeX,sizeX);
eps_temp = 1e-4;
for m=1:sizeX
for n=1:sizeX
	xpp = vecX_temp; xpp(m) += eps_temp; xpp(n) += eps_temp;
	xmp = vecX_temp; xmp(m) -= eps_temp; xmp(n) += eps_temp;
	xpm = vecX_temp; xpm(m) += eps_temp; xpm(n) -= eps_temp;
	xmm = vecX_temp; xmm(m) -= eps_temp; xmm(n) -= eps_temp;
	matH_temp(m,n) = ( funchOm_temp(xpp) - funchOm_temp(xpm) - funchOm_temp(xmp) + funchOm_temp(xmm) ) / (4.0*eps_temp*eps_temp);
end
end
echo_matH_temp = matH_temp
rcond(matH_temp)
