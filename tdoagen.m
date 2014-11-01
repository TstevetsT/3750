function [tdoa,zzc1,zzc2] = tdoagen(sigm,k,z,zc0,D,vz,dt)
%  Generates k component vector of tdoas
% sigm - Standard deviation of time measurement
% k - number of measurements
% z - True position of emitter
% D - Distance separating collectors
% vz - Vector collector velocity
% dt - Time between TDOA measurements
% tdoa - k-component vector of simulated TDOA measurements
% Written by HH Loomis on 7/7/99
% Based on nrgeogrf.m version 3.41

% s - half separation
% zz - Matrix of collector positions, 2xk

global c
disp(' Generating')
disp(k)
disp(' TDOA Observations.')
s=D/2;
zz=zeros(2,k);  

for j=1:k   
   
%r1=z-(zz(:,j)-[s;0])
%r2=z-(zz(:,j)+[s;0])
%TDOA=(1.c)*(norm(r1)-norm(r2)) + noise effect
			zz(:,j)=zc0+vz*(j-1)*dt;
			zzc2(:,j)=zz(:,j)-[s;0];
			zzc1(:,j)=zz(:,j)+[s;0];
          tdoa(j,1)=(1./c)*(norm(z-(zzc1(:,j))) - norm(z-(zzc2(:,j)))) + sigm*randn;   
 
end   

