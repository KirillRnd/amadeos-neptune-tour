function  [T,Y] = complexIntegration(y0, dV, tdV,tspan_init,lsp)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% options = odeset('AbsTol',1e-10);
Rrel=[0; 0; 0;];
RN=24622000;
dist=RN;
options = odeset('Events',@(t, y)stopByDist(t,y,Rrel,dist));
tspan_new=[tspan_init(1), tdV, tspan_init(end)];
T=zeros(0);
Y=zeros(0,6);
y0_new=y0;
for i = 1:length(tspan_new)-1
    t_start=tspan_new(i);
    t_end=tspan_new(i+1);
    tspan = linspace(t_start, t_end,lsp(i));
    if lsp(i)< 1000
        options = odeset(options,'AbsTol',1e-3);
    else
        options = odeset(options,'AbsTol',1e-10);
        options = odeset(options,'RelTol',1e-6);
    end
    if i == length(tspan_new)-1
        [t, y] = ode45(@(t,y) partialIntegration(t,y,'Neptune'),tspan,y0_new,options);
    else
        [t, y] = hashedIntegration(@(t,y) partialIntegration(t,y,'Neptune'),tspan,y0_new,options);
    end
    
    T = cat(1,T,t);
    Y = cat(1,Y,y);
    y0_new=y(end,:);
    %Учитываем импульс
    if i+1>1 && i+1<length(tspan_new)
        y0_new(4:6)=y0_new(4:6)+dV(i,:);
    end
end

end

