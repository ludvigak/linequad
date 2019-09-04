function [val, idx] = mink(array, n)

[arr_sort, idx_sort] = sort(array);
val = arr_sort(1:n);
idx = idx_sort(1:n);