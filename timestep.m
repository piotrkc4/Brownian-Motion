function [t] = timestep(startValue,endValue,nElements)
stepSize = (endValue-startValue)/(nElements-1);
t = [startValue:stepSize:endValue];
end