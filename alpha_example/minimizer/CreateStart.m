function [] = CreateStart(Pts, numVars, numPts, fileLocation)

% open file
ST = fopen(fileLocation, 'w');

% print the number of start points
fprintf(ST, '%d\n\n', numPts);

% loop through to print the points
for j = 1:numPts
    for k = 1:numVars
        fprintf(ST, '%.15e %.15e\n', real(Pts(j,k)), imag(Pts(j,k)));
    end;
    fprintf(ST, '\n');
end;

% close file
fclose(ST);


