function funcs = bertiniFormatFunctions(functions)
	funcs = [];
	[nrows ncols] = size(functions);
	for i=1:ncols
		if functions(1,i) ~= 0.0
			fprintf('F%d = (%s) - eps%d',i,char(vpa(functions(i))),i)
		end;
		fprintf(';\n')
	end;
return; 