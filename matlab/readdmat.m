function [ M ] = readdmat( filename )
  % Tao Du
  % Nov 7, 2014

  % Read a matrix from .DMAT file. The file format of a .DMAT file:
  % #col #row
  % matrix data.
  
  % Open file.
  f = fopen(filename);
  
  % Read column # and row #.
  nums = fscanf(f, '%d', 2);
  ncol = nums(1);
  nrow = nums(2);
  
  % Read matrix data.
  M = fscanf(f, '%f', [nrow, ncol]);
  
  % Close file.
  fclose(f);
end

