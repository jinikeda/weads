function [hNew,Thresholds,Param,Parent]=qp_plot_pnt(hNew,Parent,Param,data,Ops,Props)
%QP_PLOT_PNT Plot function of QuickPlot for point data sets.

%----- LGPL --------------------------------------------------------------------
%                                                                               
%   Copyright (C) 2011-2018 Stichting Deltares.                                     
%                                                                               
%   This library is free software; you can redistribute it and/or                
%   modify it under the terms of the GNU Lesser General Public                   
%   License as published by the Free Software Foundation version 2.1.                         
%                                                                               
%   This library is distributed in the hope that it will be useful,              
%   but WITHOUT ANY WARRANTY; without even the implied warranty of               
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
%   Lesser General Public License for more details.                              
%                                                                               
%   You should have received a copy of the GNU Lesser General Public             
%   License along with this library; if not, see <http://www.gnu.org/licenses/>. 
%                                                                               
%   contact: delft3d.support@deltares.nl                                         
%   Stichting Deltares                                                           
%   P.O. Box 177                                                                 
%   2600 MH Delft, The Netherlands                                               
%                                                                               
%   All indications and logos of, and references to, "Delft3D" and "Deltares"    
%   are registered trademarks of Stichting Deltares, and remain the property of  
%   Stichting Deltares. All rights reserved.                                     
%                                                                               
%-------------------------------------------------------------------------------
%   http://www.deltaressystems.com
%   HeadURL: https://svn.oss.deltares.nl/repos/delft3d/trunk/src/tools_lgpl/matlab/quickplot/progsrc/private/qp_plot_pnt.m 
%   Id: qp_plot_pnt.m 7992 2018-01-09 10:27:35Z mourits 

T_=1; ST_=2; M_=3; N_=4; K_=5;

FirstFrame=Param.FirstFrame;
Quant=Param.Quant;
Units=Param.Units;
if ~isempty(Units)
    PName=sprintf('%s (%s)',Quant,Units);
else
    PName=Quant;
end
TStr=Param.TStr;
Selected=Param.Selected;
multiple=Param.multiple;
NVal=Param.NVal;
stn=Param.stn;

DimFlag=Props.DimFlag;
Thresholds=[];

if isfield(data,'XYZ')
    X = data.XYZ(1,:,1)';
    Y = data.XYZ(1,:,2)';
    if size(data.XYZ,4)>2
        Z = data.XYZ(1,:,3)';
    else
        Z = [];
    end
elseif isfield(data,'XY')
    X = data.XY(:,1)';
    Y = data.XY(:,2)';
    Z = [];
else
    if isfield(data,'X')
        X = data.X;
    else
        X = [];
    end
    if isfield(data,'Y')
        Y = data.Y;
    else
        Y = [];
    end
    if isfield(data,'Z')
        Z = data.Z;
    else
        Z = [];
    end
end

switch NVal
    case 0
        if strcmp(Ops.facecolour,'none')
            if ishandle(hNew)
                set(hNew,'xdata',X, ...
                    'ydata',Y);
            else
                hNew=line(X,Y, ...
                    'parent',Parent, ...
                    Ops.LineParams{:});
                set(Parent,'layer','top')
            end
        else
            if ~FirstFrame
                delete(hNew)
            end
            vNaN=isnan(X);
            if any(vNaN)
                bs=findseries(~vNaN);
            else
                bs=[1 length(vNaN)];
            end
            for i=1:size(bs,1)
                if X(bs(i,1))==X(bs(i,2)) && Y(bs(i,1))==Y(bs(i,2))
                    % this patch should not influence color scaling.
                    % however, the default "1" cdata will do so
                    % we cannot set the cdata to [] immediately
                    % so, we change it after having set all color options
                    hNew(i)=patch(X(bs(i,1):bs(i,2)), ...
                        Y(bs(i,1):bs(i,2)), ...
                        1, ...
                        'edgecolor',Ops.colour, ...
                        'facecolor',Ops.facecolour, ...
                        'linestyle',Ops.linestyle, ...
                        'linewidth',Ops.linewidth, ...
                        'marker',Ops.marker, ...
                        'markersize',Ops.markersize, ...
                        'markeredgecolor',Ops.markercolour, ...
                        'markerfacecolor',Ops.markerfillcolour, ...
                        'cdata',[], ...
                        'parent',Parent);
                else
                    hNew(i)=line(X(bs(i,1):bs(i,2)), ...
                        Y(bs(i,1):bs(i,2)), ...
                        'parent',Parent, ...
                        Ops.LineParams{:});
                end
            end
            set(Parent,'layer','top')
        end
        if ~isempty(stn)
            PName = [PName ': ' stn];
        end
        qp_title(Parent,{PName,TStr},'quantity',Quant,'unit',Units,'time',TStr)
    case 1
        axestype = strtok(Ops.axestype);
        if strcmp(axestype,'Distance-Val') || strcmp(axestype,'X-Val') || strcmp(axestype,'Time-Val') || strcmp(axestype,'Time-Z')
        %if multiple(T_)
            switch axestype
                case {'Distance-Val','X-Val'}
                    x = X;
                    xdate = 0;
                otherwise
                    x = data.Time;
                    xdate = Props.DimFlag(T_)~=5;
            end
            if FirstFrame
                hNew=line(x,data.Val, ...
                    'parent',Parent, ...
                    Ops.LineParams{:});
                if xdate
                    tick(Parent,'x','autodate')
                end
            else
                set(hNew,'xdata',x,'ydata',data.Val);
            end
            if ~isempty(stn)
                Str=stn;
            else
                Str='';
            end
            qp_title(Parent,Str,'quantity',Quant,'unit',Units,'time',TStr)
        else
            switch Ops.presentationtype
                case 'values'
                    hNew=gentextfld(hNew,Ops,Parent,data.Val,X,Y);
                    
                case 'markers'
                    hNew=genmarkers(hNew,Ops,Parent,data.Val,X,Y);
                    
                otherwise
                    if ~FirstFrame
                        delete(hNew)
                    end
                    vNaN=isnan(data.Val);
                    if any(vNaN)
                        bs=findseries(~vNaN);
                    else
                        bs=[1 length(vNaN)];
                    end
                    fill = ~strcmp(Ops.facecolour,'none');
                    for i=1:size(bs,1)
                        from=bs(i,1);
                        to=bs(i,2);
                        ecol='flat';
                        fcol='none';
                        if fill && X(from)==X(to) && ...
                                Y(from)==Y(to)
                            ecol='none';
                            fcol='flat';
                            vl=from;
                        elseif from>1
                            from=from-1;
                            X(from)=NaN;
                            Y(from)=NaN;
                            data.Val(from)=NaN;
                            vl=from:to;
                        else
                            to=to+1;
                            X(to)=NaN;
                            Y(to)=NaN;
                            data.Val(to)=NaN;
                            vl=from:to;
                        end
                        hNew(i)=patch(X(from:to), ...
                            Y(from:to), ...
                            data.Val(vl), ...
                            'edgecolor',ecol, ...
                            'facecolor',fcol, ...
                            'linestyle',Ops.linestyle, ...
                            'linewidth',Ops.linewidth, ...
                            'marker',Ops.marker, ...
                            'markersize',Ops.markersize, ...
                            'markeredgecolor',Ops.markercolour, ...
                            'markerfacecolor',Ops.markerfillcolour, ...
                            'parent',Parent);
                    end
            end
            set(Parent,'layer','top')
            if strcmp(Ops.colourbar,'none')
                qp_title(Parent,{PName,TStr},'quantity',Quant,'unit',Units,'time',TStr)
            else
                qp_title(Parent,{TStr},'quantity',Quant,'unit',Units,'time',TStr)
            end
        end
    case {2,3}
        [hNew,Thresholds,Param,Parent]=qp_plot_default(hNew,Parent,Param,data,Ops,Props);
    case 4
        hNew=gentextfld(hNew,Ops,Parent,data.Val,X,Y);
end
