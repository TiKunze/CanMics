function data = fun_integrate( time,funcs,p_stream,n_streams,initvals )
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

        p_ext = p_stream(currstep);        
        n_py = n_streams(1,currstep);
        n_ei = n_streams(2,currstep);
        n_ii = n_streams(3,currstep);
        temp = data(:,currstep) + funcs(0,data(:,currstep),p_ext,n_py,n_ei,n_ii) * integ_stepsize;
        temp2 = funcs(0,temp,p_ext,n_py,n_ei,n_ii);

        data(:,currstep+1) = data(:,currstep) + (funcs(0,data(:,currstep),p_ext,n_py,n_ei,n_ii) + temp2)/2 * integ_stepsize;

    end

end

