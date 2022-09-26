function Aval = afun(i,j,x,dk,v0,dcorr)
   Aval = zeros(2*length(i),2*length(j));
   targets = x(:,i);
   sources = x(:,j);
   
   Aval = aker(targets,sources,dk,v0);
   
   
   [I,J] = ndgrid(i,j);
   tmp = Aval(1:2:end,1:2:end);
   tmp(I==J) = tmp(I==J)+1-v0*dcorr;   
   Aval(1:2:end,1:2:end) = tmp;
   
   tmp = Aval(2:2:end,2:2:end);
   tmp(I==J) = tmp(I==J)+1-v0*dcorr;  
   Aval(2:2:end,2:2:end) = tmp;
   
end