%  go2dtdoa.m
%	October 1996, Titled: "Geolocation of Electromagnetic Emitters"
%  Current version is GEOLOCHHL.doc
%  Adapted from nrgeogrf.m, version 3.41
%  Modified by hhl on 4/10/99 to provide graphical output with elipse
%  Modified by hhl on 10/5/99 to plot isochrons optionally for measured or
%    estimated-position TDOAs
%  modified by hhl on 4/1/01 (4.84) to print posit, ellipse paras, rms residuals on figure
%  modified by hhl on 6/19/01 (4.85) to fix bug in printing
%  modified by hhl on 10/14/04 (4.86) to add major axis azimuth to figure.
%  modified by hhl on 1/13/05 (4.9) to permit selection of sample interval
%  Functions called:
%      
%    tdoagen.m     generates k tdoa measurements
%    nrtdoa.m      Newton-Raphson geoposition estimator
%      elipa.m        (confidence elipse by Mark Olson)
%
%  Production version 5.1

%  Variable Declarations

% z - True position of emitter
% [x,y] - grid of x,y vectors for contour plotting
% vz - Vector collector velocity
% zz - Matrix of collector midpoint positions, 2xk
% D - Distance separating collectors
% s - half separation
% sigm - Standard deviation of time measurement
% k - number of measurements
% m - k-component vector of simulated TDOA measurements
% zi - iterated estimated position of emitter
% W - Conditioning matrix 
% mi - vector of TDOAs for current zi
% A - kx2 PARTIAL DERIV. MATRIX,
% P=inv(A'*W*A), the position covariance matrix, 2x2

%	Initialization Section

clear variables

global c 

version='5.1'  %tweak contour plots
c = 3.e08;		%Speed of Light, m/s   
dt = 1.0; 		%Time between TDOA measurements, sec
zc0=[0;0  ];	%Center of collector array at time zero
firsttime=1;      % Initializes figure 1 first time thru
isos='n';
width='n';
format compact
%	pi = asin(1.)   
  
%	INPUT SECTION   
run=4;
while run~=0;
   if run >= 4
      clear z;
      z=input(' Input True position, [x;y], m., R: ');
    end
   if run >= 3
      disp(' Center of collector cluster assumed to start at' );
      disp(zc0)
      vz=input(' Input Velocity of Collector, [vx; vy], m/sec. R: ');
      D=input(' Input Separation of collector, in X, m.: ');
      s=D/2;
   end
   if run>=2
      sigm=input(' Input Std Dev. of generated M(k), (in sec) R: ');
      dt=input(' Input sample interval (in sec) I: '); %4.9
      disp('system will run following number of measurements with ')
      disp(dt)
      disp('seconds between measurements')
      k=input(' Input number of TDOA Measurements, (<= 10) I: ');
 
%	generate value of m(k)   

[m,zzc1,zzc2]= tdoagen(sigm,k,z,zc0,D,vz,dt);

 disp('TDOAs')
 disp(m)

   end

%	OBTAIN INITIAL GUESS

% zi - iterated estimated position of emitter
% W - Conditioning matrix 
% mi - vector of TDOAs for current zi
% A - kx2 PARTIAL DERIV. MATRIX,  

   zi=input(' Input initial trial [x;y], m R: ');
% Call the Newton-Raphson Geolocation routine
   
[mi, zz, P] = nrtdoa(k,m,zzc1,zzc2,vz,sigm,zi);
   
% format long g                          

%covariance is P
% P														%%%removed squaring of diag elements (MO 3/99)
% debugged 4/6/99


% In text, for a  given ellipse probability =pc= 1 - Pe, the constant k is given by 
%      k = -2 ln Pe 
% This gives results which agree with the tabular values I have been using: 
%     Pe = 1%    ==> 99% confidence ellipse 
%     k = -2 ln .01 = 9.2 
%     Pe = 5%    ==> 95% pc confidence ellipse 
%     k = -2 ln .05 = 6.0 
%     k = cnt^2 

% disp('Containment Probability')
pc=.95;
kappa=-2*log(1-pc);

% Call elipse function
    cnt=sqrt(kappa);
    
    [xout, yout, smaj, smin, az] = elipa(P, cnt, zz(1), zz(2));
%     display('ellipse parameters')
    azd=az*(180/pi) ;
%     smaj
%     smin
%     azd
 
%compute and plot TDOA isochron contours
%lifted from geotdof2
%First set up extent of plot to include all features of interest

 mxlimit=min([-s zz(1,:)-smaj z(1)-smaj zzc2(1,:)-s]);
  pxlimit=1.*max([s z(1) zz(1)+smaj zzc1(1,:)+s]); % z(2) zz(2)+smaj zzc1(2,:)]);
  pylimit=1.*max([z(2) zz(2)+smaj zzc1(2,:)]);
  mylimit=-pylimit;

 range=max(pxlimit-mxlimit,pylimit-mylimit);

 [x,y]=meshgrid(mxlimit:.002*range:pxlimit,mylimit:.002*range:pylimit);

% now plot geometry
%
% ask if want to plot isochrons
isosm=input('Do you want to plot isochrons for measured TDOAs? (y or n):','s');
if(isosm(1)=='y')
   widthm=input('Do you want to plot those isochrons with +- sigma width? (y or n):','s');
end
isose=input('Do you want to plot isochrons for est posit TDOAs? (y or n):','s');
if(isose(1)=='y')
   widthe=input('Do you want to plot those isochrons with +- sigma width? (y or n):','s');
end
     
 
 if(run>=2)
    
    if(firsttime==1)
       fig=figure(1);
       firsttime=0;
    else
       fig=figure;
     end
    clf
    hold on
    if (pxlimit-mxlimit>pylimit-mylimit)
        axis([mxlimit pxlimit (-(pxlimit-mxlimit)/2) ((pxlimit-mxlimit)/2)]);
    else
        axis([(z(1)-(pylimit-mylimit)/2) (z(1)+(pylimit-mylimit)/2) mylimit pylimit]); 
    end
    %axis equal
    grid on
    xlabel('x position in meters')
    ylabel('y position in meters')
%   nfig=int2str(fig);
    tfig=input('Enter Figure Title: ','s');
    title(tfig);
    plot(z(1),z(2),'bx',zzc2(1,:),zzc2(2,:),'b*',zzc1(1,:),zzc1(2,:),'b*')
 end
 
%TDOA contours
if(isosm(1)=='y')
 for j=1:k
   % jth contour
  tb=(1/c)*(((x-(zzc1(1,j))).^2 + (y-(zzc1(2,j))).^2).^(.5) - ...
             ((x-(zzc2(1,j))).^2 + (y-(zzc2(2,j))).^2).^(.5));
  MT=max(max(tb));
  mt=min(min(tb));
  ttb=m(j,1);
  tlb=[ttb ttb];% tlb should be monotonically incr. vector of contour values
                    % or [z z] for a single contour of level z
  %tlb=[1 ttb];  %This one worked, but didn't match help  9/20/14
  Tb=contour(mxlimit:.002*range:pxlimit,mylimit:.002*range:pylimit,tb,tlb,'r-');
  if(widthm(1)=='y')
     ttb=m(j,1)+sigm;
     tlb=[ttb ttb];
     Tb=contour(mxlimit:.002*range:pxlimit,mylimit:.002*range:pylimit,tb,tlb,'g-');
     ttb=m(j,1)-sigm;
     tlb=[ttb ttb]; 
     Tb=contour(mxlimit:.002*range:pxlimit,mylimit:.002*range:pylimit,tb,tlb,'g-');     
  end 
 end
end
if(isose(1)=='y')
 for j=1:k
   % jth contour
  tb=(1/c).*(((x-(zzc1(1,j))).^2 + (y-(zzc1(2,j))).^2).^(.5) - ...
             ((x-(zzc2(1,j))).^2 + (y-(zzc2(2,j))).^2).^(.5));
  MT=max(max(tb));
  mt=min(min(tb));
  ttb=mi(j,1);
  tlb=[1 ttb];
  Tb=contour(mxlimit:.002*range:pxlimit,mylimit:.002*range:pylimit,tb,tlb,'b-');
  if(widthe(1)=='y')
     ttb=mi(j,1)+sigm;
     tlb=[1 ttb];
     Tb=contour(mxlimit:.002*range:pxlimit,mylimit:.002*range:pylimit,tb,tlb,'c-');
     ttb=mi(j,1)-sigm;
     tlb=[1 ttb];
     Tb=contour(mxlimit:.002*range:pxlimit,mylimit:.002*range:pylimit,tb,tlb,'c-');     
  end 
 end
end
% Plot est point and elipse
plot(zz(1),zz(2),'r.',xout,yout,'b-')
tspositx=num2str(zz(1));
tsposity=num2str(zz(2));
posittxt=strcat('Est. Posit=[',tspositx,';',tsposity,']');
%Print position estimate on figure above or below and to right of [mxlimit,0].
%Text will be on opposite side of x axis from estimated position.
tsp=(-sign(zz(2)))*(pylimit-mylimit)/20;
text(mxlimit,zzc1(2,1)-tsp,posittxt);
tsmaj=num2str(smaj);
tsmin=num2str(smin);
tsazd=num2str(azd);
elpstxt = strcat('smaj = ',tsmaj,', smin = ',tsmin,', azdeg = ',tsazd);
velstr=num2str(vz(1));
velstrcat=strcat(' vx = ',velstr);
%print ellipse paras below and to right of [mxlimit,0].
text(mxlimit,zzc1(2,1)-2*tsp,elpstxt);
%print rms residuals below ellipse params.
dm=m-mi;
rms_m = norm(dm);
tsresid=num2str(rms_m);
rmstxt=strcat('rms residual tdoas = ', tsresid);
text(mxlimit,zzc1(2,1)-3*tsp,rmstxt);
 %format short g
 
disp(tfig)
disp('Input 1 to redo geolocation with same m vector and same figure, ')
run=input('2 with new std dev, 3 with new velocity, 4 from scratch, 0 to exit: ')
 if run == 0 break
 end
% end run while loop
end
