function grads = gradsCalc(functions,symbol_order)
	grads = [];
	for elm = functions
		grads = [grads gradient(elm,symbol_order)];
	end; 
	[nrows ncols] = size(grads);
	for i=1:nrows
		fprintf('G%d = l0*(%s - %ss) + ',i,char(symbol_order(i)),char(symbol_order(i)))
		for j=1:(ncols-1)
			if grads(i,j) ~= 0.0
				fprintf('l%d*(%s) + ',j,char(vpa(grads(i,j))))
			end;
		end;
		if grads(i,ncols) ~= 0.0
			fprintf('l%d*(%s)',ncols,char(vpa(grads(i,ncols))))
		end;
		fprintf(';\n')
	end;
return; 