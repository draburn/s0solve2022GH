clear;
setprngstates();

if (1)
	hthavg = zeros(3,3);
	hsigma = randn(size(hthavg));
	hsigma = hsigma'+hsigma
	eig(hsigma)
	tmax = 10000
	for t = 1:tmax
		h = randn(size(hsigma)) .* hsigma;
		hthavg += h'*h;
	endfor
	hthavg /= tmax
	return;
elseif (1)
	hthavg = zeros(2,2);
	tmax = 10000
	for t = 1:tmax
		h = randn(2,2) .* [ 1, 5; 5, 1 ];
		h(1,2) = h(2,1);
		hthavg += h'*h;
	endfor
	hthavg /= tmax
	return;
endif



x = linspace(-1,1,51);
y = linspace(-1,1,51);
[ xx, yy ] = meshgrid( x, y );
vp = [ reshape(xx,1,[]); reshape(yy,1,[]) ];
fp = zeros(1,size(vp,2));

tmax = 10000
for t = 1 : tmax
	h = randn(2,2) .* [ 1, 5; 5, 1 ];
	h(1,2) = h(2,1);
	fp += (sum( vp .* (h*vp), 1 )).^2;
endfor
fp /= tmax;
ff = reshape( fp, size(xx) );

contourf(xx,yy,(ff).^0.1)
axis equal
axis equal
