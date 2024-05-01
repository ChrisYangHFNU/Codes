function bigscreen(arg)
%BIGSCREEN  Set graphics properties for large audiences.
%   BIGSCREEN with no arguments sets several Handle Graphics font
%   and line sizes to values appropriate for big screen projectors.
%   BIGSCREEN('off') sets them back to their "factory" defaults.

%   Copyright 2013 Cleve Moler
%   Copyright 2013 The MathWorks, Inc.

if nargin == 0 || isequal(arg,'on')
   % Turn on bigscreen
   s = get(0,'screensize');
   w = 1024;
   h = 768;
   fs = 14;
   lw = 2;
   if s(3) == 1024
      w = 768;
      h = 672;
      fs = 12;
   end
   if s(3) == 800
      w = 680;
      h = 510;
      fs = 12;
   end
   x = (s(3)-w)/2;
   y = (s(4)-h)/2;
   set(0,'defaultfigureposition',[x y w h]);
   set(0,'defaultuicontrolfontsize',fs)
   set(0,'defaulttextfontsize',fs)
   set(0,'defaultaxesfontsize',fs)
   set(0,'defaultlinemarkersize',fs)
   set(0,'defaultlinelinewidth',lw)
   set(0,'defaultuicontrolfontweight','bold')
   set(0,'defaulttextfontweight','bold')
   set(0,'defaultaxesfontweight','bold')

elseif isequal(arg,'off')
   % Turn off bigscreen
   set(0,'defaultfigureposition','factory');
   set(0,'defaultuicontrolfontsize','factory');
   set(0,'defaulttextfontsize','factory');
   set(0,'defaultaxesfontsize','factory');
   set(0,'defaultlinemarkersize','factory');
   set(0,'defaultlinelinewidth','factory');
   set(0,'defaultuicontrolfontweight','factory');
   set(0,'defaulttextfontweight','factory');
   set(0,'defaultaxesfontweight','factory');

else
   error('Unfamiliar argument')
end
