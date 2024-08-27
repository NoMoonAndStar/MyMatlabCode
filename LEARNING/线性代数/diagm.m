function A = diagm(varargin)
    A = cast([], class(varargin{1}));

    for i = 1:nargin
        A1 = varargin{i}; [n, m] = size(A);
        [n1, m1] = size(A1); A(n + 1:n + n1, m + 1:m + m1) = A1;
    end

end
