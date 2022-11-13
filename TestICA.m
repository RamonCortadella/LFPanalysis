close('all')
t=[1:100*100]/100;
f1 = 50;
f2 = 1;
f3 = 10;
f4= 20;
f5 = 25;
f6=4;
f7= 28;
f8 = 55;

IC1 = 100*repmat(sin(2*pi*f1*t),256,1)';
IC2 = [];
IC3 = [];
for i = [1:256]
    IC2(:,i) = sin(2*pi*f2*t)*i;
    IC3(:,i) = -sin(2*pi*f3*t)*abs(i-257);
    
    IC5(:,i) = sin(2*pi*f5*t)*i;
    IC6(:,i) = -sin(2*pi*f6*t)*abs(i-257);
end
IC4 = zeros(size(IC1));
IC4(:,150) = 8000*sin(2*pi*f4*t);
IC7 = zeros(size(IC1));
IC7(:,155) = 5000*sin(2*pi*f7*t);
IC8 = zeros(size(IC1));
IC8(:,170) = 2000*sin(2*pi*f8*t);
LfpSim =  IC1 + IC2+ IC3 +IC4 +IC5 +IC6 + IC7+IC8;


NumICs=10;
[ICAsignal,A,W] = fastica(LfpSim','numOfIC',NumICs);

[bX, bY]= meshgrid([1:16]',[1:16]');
v = reshape(bX,[],1);
w = reshape(bY,[],1);
figure()
szA = size(A);
for i = 1:szA(2)
    subplot(floor(sqrt(szA(2)))+1,floor(sqrt(szA(2)))+1,i)  %dim A = (nCh rows x NumICs cols)
    F=  scatteredInterpolant(v,w,A(:,i));
    bF = F(bX,bY);
    pcolor(bF)
%     caxis([-0.3 0.3])
    colorbar
end


figure()
hold on
plot(ICAsignal(4,:))
plot(ICAsignal(5,:))


