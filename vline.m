function hhh=vline(x,in1,in2,varargin)


%% Copyright notice
% By Brandon Kuczenski for Kensington Labs (brandon_kuczenski@kensingtonlabs.com) 8 November 2001

% Copyright (c) 2001, Brandon Kuczenski
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

%%

if length(x)>1  % vector input
    for I=1:length(x)
        if nargin ==1
            linetype='r:';
            label='';
        elseif nargin ==2
            if ~iscell(in1)
                in1={in1};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            label='';
        elseif nargin > 2
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            if I>length(in2)
                label=in2{end};
            else
                label=in2{I};
            end
        end
        h(I)=vline(x(I),linetype,label);
    end
else
    if nargin ==1
        linetype = 'r:';
        label    = '';
    elseif nargin == 2 
        linetype = in1;
        label    = '';
    elseif nargin > 2
        linetype = in1;
        label    = in2;
    end
    
    g=ishold(gca);
    hold on

    y=get(gca,'ylim');
    h=plot([x x],y,linetype);
    if length(label)
        xx     = get(gca,'xlim');
        xrange = xx(2)-xx(1);
        xunit  = (x-xx(1))/xrange;
        %if xunit<0.8
            text(x+0.01*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'),'rotation',90)
        %else
        %    text(x-.05*xrange, y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'),'rotation',90)
        %end
    end     

    if g==0
    hold off
    end
    if nargin==4
       handlevisibility = varargin{1};
    else
       handlevisibility ='off';
    end
    set(h,'tag','vline','handlevisibility',handlevisibility);
    % this last part is so that it doesn't show up on legends
    % NOT ALWAYS A GOOD IDEA AS THEN THE LINE DOES NOT DISAPPEAR WHEN CALLING CLA
end % else

if nargout
    hhh=h;
end
