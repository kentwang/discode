function result = if_singular(W)
% determine if W has zero rows or zero columns
	
	if min(sum(W, 2)) == 0 | min(sum(W, 1)) == 0
		result = true
	else
		result = false
	end
end