
Attempting to create a framework for density dependence in Comapdre
-------------------------------------------------------------------

This is pretty experimental, but the idea is to store matrix elements (which may be constants *or* expressions) in a long data format, then use a couple functions to take a user- or database supplied population vector and generate a density dependent matrix. Outputs are then generated via iteration.

### Density-dependent matrices

These are now implemented. `iterate_dd_mat()` can handle both user-supplied and data base matrices. Additionally, `make_mat_exprs()` is now smart enough to know when to wrap elements in calls to `eval_tidy()` and `quo()` so that end users and programmers do not need to fully understand how/why these are being used.

`CompadreDDM` matrices do not look like other matrices stored in Compadre. They are lists with 2 elements: a `data_list` which contains values for each parameter and a `mat_exprs` list which contains expressions to calculate density dependent vital rates (e.g. survival, growth, reproduction). Additionally, the `mat_exprs` list contains an expression for the matrix. See the example below for more details.
