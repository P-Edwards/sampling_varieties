function reg_4d = reg4D(points,table)
	one = points(1);
	two = points(2);
	three = points(3);
	four = points(4);
	mat = [table(one,one) table(one,two) table(one,three) table(one,four)  1;
		   table(one,two) table(two,two) table(two,three) table(two,four)  1;
		   table(one,three) table(two,three) table(three,three) table(three,four) 1;
		   table(one,four) table(two,four) table(three,four) table(four,four) 1;
		   1 1 1 1 0;];
	reg_4d = det(mat);
return; 