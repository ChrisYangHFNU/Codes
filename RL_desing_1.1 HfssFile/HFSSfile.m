%% HFSS file generation (tested with HFSS v15)
% This script modifies the coordinates of already existing HFSS file
% (TEMPLATE) containing an arbitrary planar polygon with N vertices.
%
% Input ASCII file consist of N coordinates represented by X and Y vectors
% defining the contour of the planar Rotman Lens (RL).
% To get input files run _polygon.m_ script first.

%% ACKNOWLEDGEMENT
% If you've found this code useful for your work, indicate the authors in
% your publications please:
%%
%  Dr.Michal Pokorny [pokorny.michal@seznam.cz]
%  Prof. Zbynek Raida [raida@feec.vutbr.cz]

%% 
%  SIX Research Center
%  Brno University of Technology
%  Technická 12
%  CZ-616 00 Brno
%  Czech Republic

%% Init
close all;
clear all;
clc
%% Import Rotman Lens data
DATA = dlmread('RL_parameters.tab');
N=DATA(3);  %get number of vertices
Nb=DATA(4); %get number of beam ports. It is needed due to symmetry issues.

DATA = dlmread('RL_XY_coordinates_in_mm.tab');
X=DATA(:,1);
Y=DATA(:,2);

Z=zeros(N,1);


%% Open input and ouput files
tempF = fopen('RL_template_165v.hfss','r'); %specify existing TEMPLATE file name, polygon of N vertices must be already included.
modelF = fopen('RL.hfss','w');     %specify ouput file name

%% Files processing
%Copy lines of template file until the line, where the HFSS variables could
%be written, is found
while 1                                      
       line=fgetl(tempF);   
       if strcmp(line,'			$end ''DesignDatasets'''),   break,   end
       fprintf(modelF,[line,'\n']);
end;
fprintf(modelF,[line,'\n']);

%Write polygon coordinates as HFSS variables
line='$begin ''Properties''\n';
fprintf(modelF,line);
% line=['VariableProp(''th_s'',''UD'', '''', ''',num2str(th_s*1000),'mm'')\n'];
% fprintf(modelF,line);
% line=['VariableProp(''th_cu'',''UD'', '''', ''',num2str(th_cu*1000),'mm'')\n'];
% fprintf(modelF,line);

for i=1:N
lineX = ['VariableProp(''X',num2str(i),''',''UD'', '''', ''',num2str(X(i)),'mm'')\n'];
lineY = ['VariableProp(''Y',num2str(i),''',''UD'', '''', ''',num2str(Y(i)),'mm'')\n'];
lineZ = ['VariableProp(''Z',num2str(i),''',''UD'', '''', ''',num2str(Z(i)),'mm'')\n'];
fprintf(modelF,lineX);
fprintf(modelF,lineY);
fprintf(modelF,lineZ);
end

line='$end ''Properties''\n';
fprintf(modelF,line);

%If any variables were already presented, they will be skipped
while 1
       line=fgetl(tempF);
       if strcmp(line,'			SnapMode=31'),
           fprintf(modelF,[line,'\n']); break;
       end;           
end;

%Copy lines of template file until the line, where the polygon definition
%could be written, is found
while 1
       line=fgetl(tempF);
       if strcmp(line,'										$begin ''PolylinePoints'''),   break,   end
            fprintf(modelF,[line,'\n']);
end;
fprintf(modelF,[line,'\n']);

%Write polygon definition using HFSS variables
for i=1:N
 line=['$begin ''PLPoint''\nX=''X',num2str(i),'''\nY=''Y',num2str(i),'''\nZ=''Z',num2str(i),'''\n$end ''PLPoint''\n'];
 fprintf(modelF,line);
end

for i=N-1*~mod(Nb,2):-1:1
 line=['$begin ''PLPoint''\nX=''X',num2str(i),'''\nY=''-Y',num2str(i),'''\nZ=''Z',num2str(i),'''\n$end ''PLPoint''\n'];
 fprintf(modelF,line);

end
    %    repeat writing of the first coordinate to close the polygon
% line=['$begin ''PLPoint''\nX=''X',num2str(1),'''\nY=''-Y',num2str(1),'''\nZ=''Z',num2str(1),'''\n$end ''PLPoint''\n'];
% fprintf(modelF,line);

%Skip already presented polygon definition
while 1
       line=fgetl(tempF);
       if strcmp(line,'										$end ''PolylinePoints'''),
           fprintf(modelF,[line,'\n']); break;
       end;           
end;

%Copy rest of the template, line by line
while 1
    line=fgets(tempF);
    if ~ischar(line),   break,   end
   
    j=1; %the backslash letter copier
    line2='';
     for i=1:1:length(line)
         line2(j)=line(i);
         if strcmp(line(i),'\')
             j=j+1;
             line2(j)='\';
         end;
         j=j+1;
     end;
     
    fprintf(modelF,line2);
end;

%% Close files
fclose(tempF);
fclose(modelF);
