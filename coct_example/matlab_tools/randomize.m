function output_functions = randomize(functions,target_number_of_functions)
	original_number_of_functions = size(functions);
	original_number_of_functions = original_number_of_functions(2);
	rands = rand(target_number_of_functions,original_number_of_functions);
	output_functions = [];
	for i=1:target_number_of_functions
		output_functions = [output_functions dot(functions,rands(i,:))];
	end;
return; 
