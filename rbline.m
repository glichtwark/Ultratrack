function [p1,p2]=rbline(varargin)
%function to draw a rubberband line and return the start and end points
%Usage: [p1,p2]=rbline;     uses current axes
% or    [p1,p2]=rbline(h);  uses axes refered to by handle, h

% Created by Sandra Martinka March 2002
% Edited by Dominic Farris 06/04/2016

switch nargin
case 0
  h=gca;
case 1
  h=varargin{1};
  axes(h);
otherwise
  disp('Too many input arguments.');
end
set(gcf,'Pointer','crosshair')
cudata=get(gcf,'UserData'); %current UserData
hold on;
k=waitforbuttonpress;
p1=get(h,'CurrentPoint');       %get starting point
p1=p1(1,1:2);                   %extract x and y
lh=plot(p1(1),p1(2),'r+:','LineWidth',2);      %plot starting point
udata.p1=p1;
udata.h=h;
udata.lh=lh;
set(gcf,'UserData',udata,'WindowButtonMotionFcn',{@wbmf,udata},'DoubleBuffer','on');
k=waitforbuttonpress;
p2=get(h,'Currentpoint');       %get end point
p2=p2(1,1:2);                   %extract x and y
set(gcf,'UserData',cudata,'WindowButtonMotionFcn','','DoubleBuffer','off'); %reset UserData, etc..
delete(lh);
set(gcf,'Pointer','arrow')
end

function wbmf(hObject,callbackdata,utemp)
%window motion callback function

%utemp=get(gcf,'UserData');
ptemp=get(utemp.h,'CurrentPoint');
ptemp=ptemp(1,1:2);
set(utemp.lh,'XData',[utemp.p1(1),ptemp(1)],'YData',[utemp.p1(2),ptemp(2)]);
end