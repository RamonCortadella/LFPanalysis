A = out.A;
nRows = T.nRows(ind(iFile));
nCols = T.nCols(ind(iFile));
display(size(A),'sizeA ******')
NumICs = 15;
[bX, bY]= meshgrid([1:nCols]',[1:nRows]');
v = reshape(bX,[],1);
w = reshape(bY,[],1);
figure()
for i = 1:NumICs

    display(NumICs,'NumICs ******')
    subplot(floor(sqrt(NumICs))+1,floor(sqrt(NumICs))+1,i)
    F=  scatteredInterpolant(v,w,A(:,i));
    bF = F(bX,bY);
    pcolor(bF)
    caxis([-2 2])
    colorbar()
end