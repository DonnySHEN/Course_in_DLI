close all

r = linspace(0,5,201);
theta = linspace(0,2*pi,201);
[R,T] = meshgrid(r,theta);
X = R.*cos(T);
Y = R.*sin(T);
f = @(x,y)  x.^2 + x.*y + y.^2;
fx = @(x,y) (2*x+y)/(2*sqrt(4*x.^2+4*x.*y+y.^2));
fy = @(x,y) (2*y+x)/(2*sqrt(4*y.^2+4*x.*y+x.^2));
surf(X,Y,f(X,Y)), shading interp

figure, contourf(X,Y,f(X,Y),20), colorbar, axis equal, hold on
%quiver(X,Y,2*X+Y,2*Y+X)
quiver(1,1,-fx(1,1),-fy(1,1),'r','linewidth',2)


r = linspace(0,1,201);
theta = linspace(0,2*pi,201);
[R,T] = meshgrid(r,theta);
X = R.*cos(T);
Y = R.*sin(T);
f = @(x,y)  exp(x.^2 + x.*y + y.^2);
fx = @(x,y) 2*x+y;
fy = @(x,y) 2*y+x;
figure, surf(X,Y,f(X,Y)), shading interp

figure, contour(X,Y,f(X,Y)), grid on, axis equal