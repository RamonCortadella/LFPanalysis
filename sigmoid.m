function Y = sigmoid(X,TimeScaleGlue)

Y = 1./(1+exp(-X/TimeScaleGlue));


