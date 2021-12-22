function f2 = foo0514_vizcrit__f2( x, y, pts )
	f2 = 100.0 + 0.1*y.^2;
	%f2 = y;
	%f2 = 100 + 0.1*( x.^2 + y.^2 );
	for n=1:size(pts,1)
		f2 += pts(n,5)*exp(-( (x-pts(n,1)).^2 + (y-pts(n,2)).^2 )/( 2.0*(pts(n,3).^2) ));
	end
return;
end
