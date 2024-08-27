m = [50 500 5000]

for i = m
    nnz(randi([0 1], [1 i])) / i
end
