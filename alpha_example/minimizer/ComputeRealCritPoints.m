function [RealCritPts] = ComputeRealCritPoints(Pt,NumFuncs,NumVars,InputFile,StartPoints,BertiniExecutable,RealTol)
  % Pt - column vector of the point measuring distance from
  % NumFuncs - number of functions in original system
  % NumVars - number of variables in original system
  % InputFile - input file for homotopy
  % StartPoints - start points for homotopy
  % BertiniExecutable - name of bertini file
  % RealTol - tolerance for deciding reality

  % setup final_parameters -- top parameters move to Pt and bottom parameters (perturbations) move to 0
  CreateStart([Pt;zeros(NumFuncs,1)], NumFuncs+NumVars, 1, 'final_parameters');

  % run Bertini
  eval(['!',BertiniExecutable,' ',InputFile,' ',StartPoints,' > /dev/null']); % ignore output to screen from Bertini

  % read in solutions
  [n,P] = ReadInStart(NumFuncs+NumVars+1, 'finite_solutions'); % +1 since homogeneous

  % remove the bottom coordinates and find the real solutions
  P = P(:,1:NumVars);
  ImagNorms = sqrt(sum((imag(P)').^2));
  RealCritPts = real(P(ImagNorms < RealTol,:));

return;

