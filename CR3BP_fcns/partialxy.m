function omxy = partialxy(z, mu)
 
   dn1 = ((z(1) + mu)^2 + z(2)^2 + z(3)^2)^(1.5);
   dn2 = ((z(1) + mu -1)^2 + z(2)^2 + z(3)^2)^(1.5);
   tmp = z(1) - (1-mu)*(z(1)+mu) / dn1 ;
   tmp = tmp - mu*(z(1)+mu-1) / dn2;
   omxy(1) = tmp;
   tmp = z(2) - (1-mu)*z(2) / dn1 ;
   tmp = tmp - mu*z(2) /dn2;
   omxy(2) = tmp;
   tmp = -(1-mu)*z(3) / dn1 ;
   tmp = tmp - mu*z(3) /dn2;
   omxy(3) = tmp;
end
