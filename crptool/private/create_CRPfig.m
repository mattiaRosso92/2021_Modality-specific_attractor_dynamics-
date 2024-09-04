function [h_axes,h_fig]=create_CRPfig(h,xshuttle,yshuttle)
% create_CRPfig   Creates the main figure for the CRP Toolbox
%    Used by CRP Toolbox

% Copyright (c) 2008-
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1998-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2016/09/23 12:53:42 $
% $Revision: 4.12 $
%
% $Log: create_CRPfig.m,v $
% Revision 4.12  2016/09/23 12:53:42  marwan
% link axis of RP and data for simultaneous effects when zooming etc.
%
% Revision 4.11  2016/03/03 14:57:22  marwan
% updated input/output check (because Mathworks is removing downwards compatibility)
% bug in crqad_big fixed (calculation of xcf).
%
% Revision 4.10  2014/10/06 07:16:44  marwan
% bug in colorbar scaling fixed
%
% Revision 4.9  2009/03/24 08:36:04  marwan
% copyright address updated
%
% Revision 4.8  2008/02/25 11:45:27  marwan
% fix of the colorbar bug
%
% Revision 4.7  2004/11/10 07:04:28  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

global props

  scr=get(0,'ScreenSize'); 
  x=xshuttle(:,2:end);
  xscale=xshuttle(:,1);
  y=yshuttle(:,2:end);
  yscale=yshuttle(:,1);

  h_fig=figure('Tag','CRPFig',...		% Plot Figure
            'Position',[scr(3)/4 scr(4)/8 3*scr(3)/8 3*scr(4)/4 ],...
	    'Color',props.window.Color,...
            'NumberTitle','off',...
	    'Name',['Cross Recurrence Plot (' h ')'],...
	    'DeleteFcn','crp smartclose',...
	    'PaperPosition',[0.25 0.25 7.7677 11.193],...
	    'PaperType','a4',...
	    'PaperOrientation','portrait');
  set(h_fig,props.window,'Units','Norm')
	    
  % create background axis object if x and y are different
  if ~isequal(x,y)	    
     h1=axes(props.axes,'Units','norm','Position',[.1 .78 .8 .15]);	% Data1 Plot

     if size(y,2)>1
       plot(y,'--','Tag','Data1')
     else
       plot(yscale,y,'r')
     end
     set(h1,'Tag','DataPlot1',...
               'Xcolor','r','ycolor','r',...
	       'XaxisLocation','top',...
	       'YaxisLocation','right',...
	       'UserData',yshuttle)
  end
  h2=axes(props.axes,'Units','norm',...	  % Data2 Plot
           'Position',[.1 .78 .8 .15],...
	    'Xcolor','k','ycolor','k',...
	    'XaxisLocation','bottom','yaxislocation','left');

  if size(y,2)>1
    plot(x,'-','Tag','Data2')
  else
    plot(xscale,x,'k')
  end
  set(h2,'Tag','DataPlot2','UserData',xshuttle)
  if isequal(x,y)
      h1 = h2; 
  else
      set(h2,'color','none')
  end	    
  
  grid on

  if max(abs(x)) > max(abs(y))
      scaling=abs(get(h2,'ylim'));
      if scaling(1)>scaling(2)
         set(h2,'ylim',[-scaling(1) scaling(1)])
      else
         set(h2,'ylim',[-scaling(2) scaling(2)])
      end
      scaling=get(h2,'ylim');
      set(h1,'ylim',scaling)
  else
      scaling=abs(get(h1,'ylim'));
      if scaling(1)>scaling(2)
         set(h1,'ylim',[-scaling(1) scaling(1)])
      else
         set(h1,'ylim',[-scaling(2) scaling(2)])
       end
      scaling=get(h1,'ylim');
      set(h2,'ylim',scaling)
  end

  h1a=title('Underlying Time Series','units','normalized');
  h2a=get(h1a,'Position');
  set(h1a,'Position',[h2a(1) h2a(2)+.12 h2a(3)])  

  h_axes=axes(props.axes,'Units','norm',...	% CRP Plot
	    'Color',[1 1 1], ...
            'Tag','CRPPlot',...	
            'Position',[.1 .12 .8 .8*17/23]);
	  
        
  linkaxes([h1 h2 h_axes],'x')  
	    
  h1b=imagesc( 'Tag','CRPData','cdata',[]);
  h1b=title('','units','normalized');
  set(h1b,props.titletext)
  h2b=get(h1b,'Position');
  set(h1b,'Position',[h2b(1) h2b(2)-.03 h2b(3)])

  h1=text(mean(xlim(h2)), .5,'busy...',...			% Text Busy
            props.normaltext,...
	    'Tag','Status',...
            'Visible','off',...
	    'HorizontalAlignment','center',...
	    'VerticalAlignment','middle',...
	    'FontSize',18);

  colormap(french(256))
  h1=colorbar('horiz');				% Colorbar
  set(h1,'Visible','off','Position',[.1 .07 .8 .02],'Tag','Colorbar')
  set(get(h1,'children'),'Visible','off')
  set(findobj('Tag','CRPPlot'),'Position',[.1 .12 .8 .8*17/23])

  cm={'hsv';'hot';'gray';'french';'bone';'copper';...    % Colormap
         'pink';'flag';'lines';'colorcube';...
	 'jet';'prism';'cool';'autumn';...
	 'spring';'winter';'summer'};
  h0=uimenu('Label','Colormap','Tag','cm');
  for i=1:length(cm);
    h1=uimenu(h0,'Label',cm{i},'Checked','Off',...
           'Tag',num2str(i),...
           'Callback','crp colormap');
    if i==4, set(h1,'Checked','On'), end
  end
  h1=uimenu(h0,'Label','b/w','Tag','18','Callback','crp colormap');
  h1=uimenu(h0,'Label','inverse','Tag','19','Separator','On','Callback','crp colormap');

  h1=uimenu('Label','SmartClose',...           % SmartClose
         'Callback','crp smartclose');

  clear h1 h2 xshuttle yshuttle
  set(h_fig, 'HandleVis','CallBack')
