function [U, V, K_plus, L_plus] = sampleDIBP(a, b, I, J)
% Sample features from dependent Indian Buffet Process
% Initialization essentially use two individual IBP !
	[U, K_plus] = sampleIBP(a, I);
	[V, L_plus] = sampleIBP(b, J);    
end