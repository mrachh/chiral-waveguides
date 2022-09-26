function G =  get_gfunself(s0,delt,de)

G = complex(zeros(2,2));
isgn = 1;
f1 = @(x,y) fgaussint(x,y,s0,de,delt,isgn);
isgn = 0;
f2 = @(x,y) fgaussint(x,y,s0,de,delt,isgn);

r0 = 6*sqrt(2)/s0;
G(1,1) = integral2(f1,-r0,r0,-r0,r0,'Method','Iterated','AbsTol',1e-6,'RelTol',1e-7);
G(2,2) = integral2(f2,-r0,r0,-r0,r0,'Method','Iterated','AbsTol',1e-6,'RelTol',1e-7);
end