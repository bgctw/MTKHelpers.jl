```@meta
CurrentModule = MTKHelpers
```

# Reading and Writing DataFrames that contain ComponentArrays

`MTKHelpers` supports setting parameters that are symbolic arrays,
by specifying them as `ComponentVector`.

Therefore, it is often useful to store parameters as `ComponentVector` in 
one columns of a `DataFrame`. A `DataFrame` with such a columns, however, 
is not easily written and retrieved from a interoperable file, such as CSV.

The following functions write and read data of such a `DataFrame` in CSV
format that can be read by other tools. By using additional information
on the Axes stored in comments, the `ComponentVector`s can be reconstructed.

```@docs
write_csv_cv
```







