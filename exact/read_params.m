%  Reads parameters from a text file
%
%  Input variable: paramfile = string containing name of parameter file
%			       (e.g junk.par)
%  File syntax:
%    fontsize = 20  % this is numeric parameter
%    holetype = 'tunnel'  % this is a string parameter
%    % This is a comment line
%
%  WARNINGS:
%    !! no ; after the values please !!
%
%    !! sometimes the next line after a comment is not read !!
%    example:  xmin = 100  % interval for measurements
%	       xmax = 200  
%         the parameter xmax may not be read in and the program will
%	  crash with the message "unknown variable or function xmax"
%    The remedy is to leave an empty line after each line that contains
%    a comment

[ parnames, parvals ] = textread( paramfile, '%s = %s %*[^\n]', 'commentstyle', 'matlab' );

for ipar = 1:length( parnames );
    cmdeval = strcat(parnames( ipar ), ' = ', num2str(char(parvals(ipar))), ';');
    eval( char(cmdeval ));
end


