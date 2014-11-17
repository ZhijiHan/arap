function [  ] = writeoff( V, F, filename )
  % Tao Du
  % Nov 17, 2014
  % This funciton writes .off files. It takes the vertex matrix and face
  % matrix, a .off filename, then writes V and F into file.
  % 
  % The definition of vertex matrix V: it is a # by 3 matrix, each row
  % represents a vertex in the object.
  % 
  % The definition of face matrix F: it is a # by 3 matrix, each row
  % represents a triangle, the 3 numbers in the triangle are the row number
  % of the vertices in V.
  
  % Open text file.
  f = fopen(filename, 'w');
  
  % Write the header: it should be OFF.
  fprintf(f, '%s\n', 'OFF');
  
  % Vertex number, face number and edge number.
  vnum = size(V, 1);
  fnum = size(F, 1);
  enum = fnum * 3 / 2;
  fprintf(f, '%d %d %d\n', vnum, fnum, enum);
  
  % Write vertex matrix.
  fprintf(f, '%f %f %f\n', V');
  
  % Write face matrix. Decrement F by 1 because indices in matlab start
  % from 1.
  F = F - 1;
  
  % Augment a new row in F: this column contains the number of
  % vertex in each row.
  F = [3 * ones(fnum, 1) F];
  fprintf(f, '%d %d %d %d\n', F');
  
  % Close the file.
  fclose(f);
end
