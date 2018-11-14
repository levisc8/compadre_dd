## Attempting to create a framework for density dependence in Comapdre

This is pretty experimental, but the idea is to store matrix elements (which may
be constants _or_ expressions) in a long data format, then use a couple functions
to take a user- or database supplied population vector and generate a density
dependent matrix. Outputs are then generated via iteration.
