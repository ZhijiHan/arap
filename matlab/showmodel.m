function [  ] = showmodel( V, F, varargin )
  % Tao Du
  % Nov 17, 2014

  % Given vertices, faces and possibly constraints, show the model.
  mkdir display;
  writeoff(V, F, 'display/model.off');
  command = '/home/taodu/research/arap/build/demo_bin display/model.off';
  
  % If we have a third argument, extract it as the constraints.
  if ~isempty(varargin)
    C = varargin{1};
    writedmat(C, 'display/model.dmat');
    command = [command ' display/model.dmat'];
  end

  % Run C++ program to display it. 
  system(command);
  
  rmdir('display', 's');
end

