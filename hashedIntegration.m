function [t,y] = hashedIntegration(fun,tspan,y0,options)
%hashedIntegration Функция для хэширования и сохранения результатов интегрирования
%   На диск в папку hashedInt сохраняется результат интегрирования для
%   переиспользования в будущем.
initIntParams = struct;
initIntParams.tspan = tspan;
initIntParams.y0 = y0;
initIntParams.options = options;

hashInt = DataHash(initIntParams);
fileHash = strcat('hashedInt/',hashInt,'.mat');
try
    sInt = load(fileHash);
    t = sInt.t;
    y = sInt.y;
catch ME
    if (strcmp(ME.identifier,'MATLAB:load:couldNotReadFile'))
        [t, y] = ode45(fun,initIntParams.tspan,initIntParams.y0,initIntParams.options);
        save(fileHash,'t','y')
    else
        rethrow(ME);
    end   
end
end

