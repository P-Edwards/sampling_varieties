function star_d = starD(points,table)
	one = points(1);
	two = points(2);
	three = points(3);
	four = points(4);
	five = points(5);
	six = points(6); 
	mat = [table(one,one) table(one,two) table(one,three) table(one,four) table(one,five) 1;
		   table(one,two) table(two,two) table(two,three) table(two,four) table(two,five) 1;
		   table(one,three) table(two,three) table(three,three) table(three,four) table(three,five) 1;
		   table(one,four) table(two,four) table(three,four) table(four,four) table(four,five) 1;
		   table(one,six) table(two,six) table(three,six) table(four,six) table(five,six) 1;
		   1 1 1 1 1 0;];
	star_d = det(mat);
return; 