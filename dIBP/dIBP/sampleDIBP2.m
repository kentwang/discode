function [U, V, K_plus, L_plus] = sampleDIBP2(a, b, I, J)
% Sample features from dependent Indian Buffet Process
% with equal K_plus and L_plus
	K_plus = 0;
	L_plus = 1;
	while(K_plus != L_plus || K_plus <= 2)
		[U, K_plus] = sampleIBP(a, I);
		[V, L_plus] = sampleIBP(b, J);    
	end
end