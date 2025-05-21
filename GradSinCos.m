function GradSinCos=GradSinCos(m,n,x,y)
    %SinCos=sin(2*pi*m*x).*cos(2*pi*n*y);
    xderiv=-cos(2*pi*m*x)*2*pi*m.*cos(2*pi*n*y);
    yderiv=sin(2*pi*m*x).*sin(2*pi*n*y)*2*pi*n;
    GradSinCos=[xderiv;yderiv];
end