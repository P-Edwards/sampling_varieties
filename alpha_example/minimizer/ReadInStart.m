function [numPts, Pts] = ReadInStart(numVars, fileLocation)

% open file
RF = fopen(fileLocation, 'r');

% make sure the file exists
if RF == -1
    % file does not exist
    disp(['The file ''', fileLocation, ''' does not exist!']);
    numPts = 0;
    Pts = 0;
    return;
end;

% read in the number of points
numPts = fscanf(RF, '%d', 1);

% setup Pts
Pts = zeros(numPts, numVars);

% loop through to read in the points
tempR = 0; tempI = 0;
for k = 1:numPts
    A = fscanf(RF,'%e',2*numVars);
    if feof(RF)
      disp(['Ran out of points!']);
      fclose(RF);
      numPts = 0;
      Pts = 0;
      return;
    end;
    Pts(k,:) = A(1:2:end) + 1i*A(2:2:end);
end;

% close file
fclose(RF);

RF = 0;

