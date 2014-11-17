function [  ] = writedmat( M, filename )
  % Tao Du
  % Nov 17, 2014

  % Read a matrix from .DMAT file. The file format of a .DMAT file:
  % #col #row
  % matrix data.
  
  % Open file.
  f = fopen(filename, 'w');
  
  % Write column # and row #.
  ncol = size(M, 2);
  nrow = size(M, 1);
  fprintf(f, '%d %d\n', ncol, nrow);
  
  % Write matrix data.
  fprintf(f, '%f\n', M);
  
  % Close file.
  fclose(f);
end

