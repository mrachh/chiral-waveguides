function f = fgaussint(x,y,s0,de,delt,isgn)
    rr = -(x.^2 + y.^2).*s0^2/2;
    rden = (x-de).^2 -1 -1j*delt -y.^2;
    rnum = (x-de) + (-1)^isgn;
    f = exp(rr).*rnum./rden;
end