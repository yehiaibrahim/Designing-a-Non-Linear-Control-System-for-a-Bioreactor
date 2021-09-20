a= [-3:.1:3];
b= [-3:.1:3];
[x,y]=meshgrid(a,b);
z= (4*y.*(x.^2))./(0.06 + y + 0.3*(y.^2)) - 2*(x.^2) + (-5.714*x.*(y.^2))./(0.06 + y + 0.3*(y.^2))+4-(2*(y.^2))
surf(x,y,z);
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');