A=[4 -1 0 -1 0 0;-1 4 -1 0 -1 0;0 -1 4 0 0 -1;-1 0 0 4 -1 0;0 -1 0 -1 4 -1;0 0 -1 0 -1 4];
disp(A);
b=9*[7; 3; -6; 7; 3; -6 ];
disp(b);
u=linsolve(A,b);
disp(u);