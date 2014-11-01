function [mi, zz, P] = nrtdoa(k,m,zzc1,zzc2,vz,sigm,zi)
%function computes position estimate, zi, and covariance matrix, P
% vz - Vector collector velocity
% zzc1 - Matrix of collector 1 positions, 2xk
% zzc2 - Matrix of collector 2 positions, 2xk
% sigm - Standard deviation of time measurement
% zi - iterated estimated position of emitter
% Written by HH Loomis on 7/7/99
% Version 1.6

% W - Conditioning matrix 
% mi - vector of TDOAs for current zi
% A - kx2 PARTIAL DERIV. MATRIX,

global c 

%	calculate Weighting Matrix, W

 %    limit=1.4*max(z);
 %     [x,y]=meshgrid(-limit:.02*limit:limit,-limit:.02*limit:limit);
 
 W=zeros(k,k);

 for jj=1:k
	      W(jj,jj)=1/(sigm^2);
 end
% disp('W=' )
% disp(W)

%	start estimation loop with ii = 0

ii=0;
while ii<10
%  disp('Begin Iteration #:')
%  disp(ii)


%calculate mi as function of zi, initial trial value of position
  
  for jj=1:k    

    %  zz(:,jj)=zc0+vz*(jj-1)*dt;
    mi(jj,1)=(1./c)*(norm(zi-(zzc1(:,jj))) - norm(zi-(zzc2(:,jj))));
%    mmi(jj,1)=(1/c)*(sqrt((zi-(zzc1(:,jj)))'*(zi-(zzc1(:,jj)))) - sqrt((zi-(zzc2(:,jj)))'* ...
%       (zi-(zzc2(:,jj)))));
  end  
  
%	Compute rms error in mi
	 
  dm=m-mi;
  rms_m = norm(dm);
  
%   disp('mi, dm, rms_m')
%   disp(mi)
%   disp(dm)
%   disp(rms_m)
%   disp('mmi')
%   disp(mmi)
  
  if rms_m <= .1*sigm break 
  else
    
%	CALCULATE PARTIAL DERIV. MATRIX, A(K,2)
    for jj=1:k
      A(jj,:)=(1/c)*(((zi-(zzc1(:,jj)))'/norm(zi ...
         -(zzc1(:,jj)))) - ...
       ((zi-(zzc2(:,jj)))'/norm(zi-(zzc2(:,jj)))));
    end
%    disp('Matrix A')
%    disp(A)
%disp ('Condition of At*W*A')
%cond(A'*W*A)

%	COMPUTE increment in z, dzi

    P=inv(A'*W*A);
    dzi=(P)*(A'*W)*dm;

% disp('dzi')
% disp(dzi)

%	COMPUTE next estimate for zi

    zip1=zi+dzi;
% disp('zip1')
% disp(zip1)

%ii used as loop index since i is sqrt(-1) in matlab

%    disp('end of iteration loop #')
%    disp(ii)
    ii=ii+1;
    zi=zip1;
    
  end %for if -- else
  
end %while ii
zz=zi;
%end iteration loop

% disp('ii')
% disp(ii)
% disp('zi, dzi')
% disp(zi)
% disp(dzi)
% disp('dm')
% disp(dm)
disp('rms_m')
disp(rms_m)
disp ('P')
disp (P)


