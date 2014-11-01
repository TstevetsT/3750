function [xout, yout, smaj, smin, theta] = elipa(PK, c, xt, yt)

%calculates error ellipsoids given error covariance
%and estimate position
%	PK is covariance matrix
%	c is conficence region (c=3 for 98%)
%	xt is ellipse center x coord
%	yt is ellipse center y coord

%	xout is a vector containing x coord points of the ellipse
%	yout is a vector containing y coord points of the ellipse
%	smaj is the length of the semi-major axis
%	smin is the length of the semi-minor axis
%	theta is the orientation of the semi-major axis

%  to plot the ellipse, use:  plot(xout,yout)

%adapted from Stephen L. Spehn's Errellip.m (15 Nov 89)  by Mark Olson
%minor mods by HH Loomis on 7/14/99

%get eigenvalues (lam) and eigenvectors(V)
[V,lam] = eig(PK);
sigx = sqrt(lam(1,1));
sigy = sqrt(lam(2,2));

%parameterized ellipse
t = 0:2*pi/100:2*pi;
x = sigx*c*cos(t);
y = sigy*c*sin(t);

%translate to eigenvectors space and center at tgt posit
xout = x*V(1,1) + y*V(1,2) + xt;
yout = x*V(2,1) + y*V(2,2) + yt;

%report semimajor and semiminor axes lengths and orientation
% corrected by HHL on 7/14/99 by removing factor of 2 in following 2 equations.
smaj = c*max([sigx sigy]);
smin = c*min([sigx sigy]);

%theta = atan2(V(2,2),V(1,2)) + pi/2*(sigy > sigx);	%sigy > sigx ==> y is smaj
theta = atan2(V(2,1),V(1,1)) + pi/2*(sigy > sigx);	%sigy > sigx ==> y is smaj
