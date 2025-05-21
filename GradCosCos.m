function GradCosCos=GradCosCos(m,n,x,y)
    %CosCos=cos(2*pi*m*x).*cos(2*pi*n*y);
    xderiv=sin(2*pi*m*x)*2*pi*m.*cos(2*pi*n*y);
    yderiv=cos(2*pi*m*x).*sin(2*pi*n*y)*2*pi*n;
    GradCosCos=[xderiv;yderiv];
end