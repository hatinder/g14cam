
fun=@(x,y) sin(x.^2+y.^2);
resultE=integral2(fun,0,1,0,1);
disp("E : ");
disp(resultE);
ymax=@(x) x;
resultF=integral2(fun,0,1,0,ymax);
disp("F : ");
disp(resultF);

