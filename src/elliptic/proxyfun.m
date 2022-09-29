function [Kpxy,nbr] = proxyfun(xflam,slf,nbr,l,ctr,x,dk,v0,opdims,pr,pin)

   pxy = pr.*l + ctr;
   
   
   slfpts = idivide(int64(slf(:)-1),int64(opdims(2)))+1;
   [slfuni,~,islfuni] = unique(slfpts);
   islfuni2 = (islfuni-1)*opdims(2) + mod(slf(:)-1,opdims(2))+1;
   xuse = x(:,slfuni);
   [~,m] = size(pxy);
   [~,n] = size(xuse);
   Kpxy1 = zeros(2*m,2*n);
   Kpxy2 = zeros(2*m,2*n);
   Kpxy3 = zeros(2*m,2*n);
   Kpxy4 = zeros(2*m,2*n);
   dx = bsxfun(@minus,pxy(1,:)',xuse(1,:));
   dy = bsxfun(@minus,pxy(2,:)',xuse(2,:));
   dr = sqrt(dx.^2 + dy.^2);
   h0 = besselh(0,dk*dr)/4;
   Kpxy1(1:2:end,1:2:end) = h0;
   Kpxy2(2:2:end,1:2:end) = h0;
   Kpxy3(1:2:end,2:2:end) = h0;
   Kpxy4(2:2:end,2:2:end) = h0;
   Kpxy1 = Kpxy1(:,islfuni2);
   Kpxy2 = Kpxy2(:,islfuni2);
   Kpxy3 = Kpxy3(:,islfuni2);
   Kpxy4 = Kpxy4(:,islfuni2);
   
   Kpxy = [Kpxy1;Kpxy2;Kpxy3;Kpxy4];
%    
%    nbrpts = idivide(int64(nbr(:)-1),int64(opdims(1)))+1;
%    dr = (x(:,nbrpts) - ctr(:))./l;
%    nbr = nbr(sum(dr.^2) < 1.5^2);
    
end