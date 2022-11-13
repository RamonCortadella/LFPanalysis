
x = [1:100]/10+randn([1,100]);
y = [1:100]/10+randn([1,100]);
z = [1:100]/10+randn([1,100]);

figure()
MarkerSize= 100;
c = 1:length(x);
scatter3(x,y,z,MarkerSize,c,'filled')