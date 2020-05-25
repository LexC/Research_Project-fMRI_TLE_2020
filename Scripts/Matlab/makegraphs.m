function [ G ] = makegraphs( map )
% Return graphs and it's parameters
%
%   [ G ] = makegraphs( map )

u=1;
G{18}=[];

for i=0.05:0.05:0.90
    
    r=map>i;
    
    G{u}=graph(r,'OmitSelfLoops');

    u=u+1;
end

end

