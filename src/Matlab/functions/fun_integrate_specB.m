function data = fun_integrate( time,funcs,p_stream,B_stream,initvals )
%FUN_INTEGRATE Summary of this function goes here
%   Detailed explanation goes here

N=length(time);
integ_stepsize=time(2);

if length(p_stream) ~= N
    error('check p_stream length')
end

data = zeros(length(initvals),N);
data(:,1) = initvals;

    for currstep=1:N-1

        ppy = p_stream(currstep);        
        B = B_stream(currstep);
        
        temp = data(:,currstep) + funcs(0,data(:,currstep),ppy,B) * integ_stepsize;
        temp2 = funcs(0,temp,ppy,B);

        data(:,currstep+1) = data(:,currstep) + (funcs(0,data(:,currstep),ppy,B) + temp2)/2 * integ_stepsize;

    end

end

