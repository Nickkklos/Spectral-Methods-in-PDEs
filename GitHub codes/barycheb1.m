function p = barycheb1(x, f, xx)
N = length(x)-1;
w = ones(N+1,1);
w(1)   = 0.5;
w(end) = 0.5;
w = w .* (-1).^(0:N)';

p = zeros(size(xx));
for k = 1:numel(xx)
    dif = xx(k) - x;
    j = find(abs(dif) < 1e-14, 1);
    if ~isempty(j)
        p(k) = f(j);
    else
        tmp = w ./ dif;
        p(k) = (tmp.'*f) / sum(tmp);
    end
end
end