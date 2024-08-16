
"""
    write_csv_cv(con::IO, df0::DataFrame)
    read_csv_cv(con::IO)::DataFrame

Write/Read a DataFrame that contains columns of type ComponentArrays   
to a CSV file.
Each ComponentArray column is replaced by length(cv) columns,
and information on the axis is generated as a comment in front of the CSV.

When re-reading the CSV, the respecive ComponentArrays are recreated 
usingn this information and are replaced for the read plain columns.

Make sure to not duplicate names. If there is a ComponentArray column
`u0` containing three Floats, the following three columns are added 
(and should not exist before): `u0_1, u0_2, u0_3`.
"""
function write_csv_cv end,
function read_csv_cv end
