function [rhs1,rhs2] = get_rhs_11(v1,v2,de,dh)
     v1dy = 0*v1;
     v1dy(2:(end-1)) = (v1(3:end)-v1(1:(end-2)))/(2*dh);
     v2dy = 0*v2;
     v2dy(2:(end-1)) = (v2(3:end)-v2(1:(end-2)))/(2*dh);
     rhs1 = -v2dy+1i*(1-de)*v1;
     rhs2 = -v1dy-1i*(1+de)*v2;
     
end

