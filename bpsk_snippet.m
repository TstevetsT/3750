
%	Produce a random sequence of {0,1}
%rand('uniform');
j = 1;
for i=1:n;
 if fix((i-1)*fb) < fix(i*fb);
  r=rand ;
  if r <0.5;
     j=-1;
   else; 
     j=1;
  end;
 end;
 p(i)=j;
end;

