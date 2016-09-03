function regular_d = regularD(points,table)
	one = points(1);
	two = points(2);
	three = points(3);
	four = points(4);
	five = points(5);
	mat = [table(one,one) table(one,two) table(one,three) table(one,four) table(one,five) 1;
		   table(one,two) table(two,two) table(two,three) table(two,four) table(two,five) 1;
		   table(one,three) table(two,three) table(three,three) table(three,four) table(three,five) 1;
		   table(one,four) table(two,four) table(three,four) table(four,four) table(four,five) 1;
		   table(one,five) table(two,five) table(three,five) table(four,five) table(five,five) 1;
		   1 1 1 1 1 0;];
	regular_d = det(mat);
return; 