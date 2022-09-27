function Aval = aker(targets,sources,dk,v0)
   [~,m] = size(targets);
   [~,n] = size(sources);
   Aval = zeros(2*m,2*n);
   dx = bsxfun(@minus,targets(1,:)',sources(1,:));
   dy = bsxfun(@minus,targets(2,:)',sources(2,:));
   dr = sqrt(dx.^2 + dy.^2);
   h0 = besselh(0,dk*dr)/4;
   h1 = besselh(1,dk*dr)/4;
   
   h1overr = h1./dr;
   h0(dr == 0) = 0;
   h1overr(dr == 0) = 0;
    
   Aval(1:2:end,1:2:end) = -v0*(1i*dk*h0 + dk*h1overr.*dx);
   Aval(1:2:end,2:2:end) = -v0*(dk*h1overr.*dy);
   Aval(2:2:end,1:2:end) = -v0*(dk*h1overr.*dy);
   Aval(2:2:end,2:2:end) = -v0*(1i*dk*h0 - dk*h1overr.*dx);

end