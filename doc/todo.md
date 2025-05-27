# Next

- Rename to blas
- Move non-blas/numa functionality to santoku.vector.integer/number (based on kvec)
- This API should operate directly on santoku.vector.number via C API

# Consider

- Can we replace dense codes with sparse set-bits format
- Heap for top_chi2 and top_mi instead of keeping full list
