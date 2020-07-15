function [parameters,fit_result] = FRAP_fit(traces)
%assume the frist column of traces is time stamp
%fit as a*exp(b*x)+c, return parameters of a, b. 

    [m,n] = size(traces);
    fit_result = zeros(m,n);
    fit_result(:,1) = traces(:,1);
    fit_result(1,2:end) = 1;
    parameters = zeros(2,n);
    
    g = fittype('a*(1-exp(b*x))');
    
    for i = 2 : n
        f = fit(traces(2:end,1)-traces(2),traces(2:end,i),g,'StartPoint',[1,-0.1]);
        parameters(:,i) = [f.a; f.b];
        fit_result(2:end,i) = f.a*(1-exp(f.b*(traces(2:end,1)-traces(2))));
    end

end