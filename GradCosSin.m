function GradCosSin=GradCosSin(m,n,x,y)
    %CosSin=cos(2*pi*m*x).*sin(2*pi*n*y);
    xderiv=sin(2*pi*m*x)*2*pi*m.*sin(2*pi*n*y);
    yderiv=-cos(2*pi*m*x).*cos(2*pi*n*y)*2*pi*n;
    GradCosSin=[xderiv;yderiv];
end