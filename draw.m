A = load("result.txt");
N = 32;
v = reshape(A(1,:),N,N);
L = 1/N:1/N:1;
[x,y] = meshgrid(L, L);
pcolor(x, y, v);
shading flat
colorbar