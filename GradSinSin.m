function GradSinSin=GradSinSin(m,n,x,y)
    %SinSin=sin(2*pi*m*x).*sin(2*pi*n*y);
    xderiv=-cos(2*pi*m*x)*2*pi*m.*sin(2*pi*n*y);
    yderiv=-sin(2*pi*m*x).*cos(2*pi*n*y)*2*pi*n;
    GradSinSin=[xderiv;yderiv];
end