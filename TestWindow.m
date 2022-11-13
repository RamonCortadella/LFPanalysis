close('all')
TimeScaleGlue = 0.1;
Fs = 651;
WindowTemp = ones([10000,1]);
Window1 = sigmoid([1:10000]/Fs-1*TimeScaleGlue,TimeScaleGlue);

Window= WindowTemp.*Window1'.*flip(Window1');

figure()
plot([1:10000]/Fs,Window)

Window = repmat(Window,10);