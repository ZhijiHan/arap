function [ V, F ] = readoff( filename )
  % Tao Du
  % Nov 7, 2014
  % This funciton reads .off files. It takes a .off filename, then returns
  % the vertex matrix and face matrix of the object.
  % 
  % The definition of vertex matrix V: it is a # by 3 matrix, each row
  % represents a vertex in the object.
  % 
  % The definition of face matrix F: it is a # by 3 matrix, each row
  % represents a triangle, the 3 numbers in the triangle are the row number
  % of the vertices in V.
  
  % Open text file.
  f = fopen(filename, 'r');
  
  % Read the header: it should be OFF.
  header = fscanf(f, '%c', 3);
  if ~strcmp(header, 'OFF')
    disp('Invalid .OFF file.');
    fclose(f);
  end
  
  % Read vertex number, face number and edge number.
  nums = fscanf(f, '%d', 3);
  vnum = nums(1);
  fnum = nums(2);
  
  % Read vertex matrix.
  V = fscanf(f, '%f', [3, vnum]);
  V = V';
  
  % Read face matrix.
  F = fscanf(f, '%f', [4, fnum]);
  F = F';
  
  % Shift F by 1 because indices in matlab strat from 1.
  F = F + 1;
  
  % Filter out the first row in F: this column contains the number of
  % vertex in each row.
  F = F(:, 2 : end);
  
  % Close the file.
  fclose(f);
end
