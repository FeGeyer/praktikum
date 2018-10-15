a = [0 0 sqrt(2) 0 sqrt(2) 0 sqrt(2) 0 0;
1 1 1 0 0 0 0 0 0;
0 0 sqrt(2) 0 0 0 0 0 0;
0 0 0 1 1 1 0 0 0;
0 sqrt(2) 0 0 0 sqrt(2) 0 0 0;
0 0 0 0 0 0 1 1 1;
sqrt(2) 0 0 0 sqrt(2) 0 0 0 sqrt(2);
0 0 1 0 0 1 0 0 1;
0 0 0 sqrt(2) 0 0 0 sqrt(2) 0;
0 1 0 0 1 0 0 1 0;
0 0 0 0 0 0 sqrt(2) 0 0;
1 0 0 1 0 0 1 0 0;];

c = inv(transpose(a) * a)
d = inv(c)
