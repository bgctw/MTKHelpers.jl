# MTKHelpers

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bgctw.github.io/MTKHelpers.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bgctw.github.io/MTKHelpers.jl/dev)
[![Build Status](https://github.com/bgctw/MTKHelpers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bgctw/MTKHelpers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/bgctw/MTKHelpers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/bgctw/MTKHelpers.jl)

Utilities that ease the development using [ModelingToolkit.jl](https://mtk.sciml.ai/dev/).
Hopefully, these will become obsolete with further development in `ModelingToolkit`.

With the introduction of MTKParameters object in ModelingToolkit 9, the 
indexing into vectors by position is not supported any more. Hence, MTKHelpers need
to be adapted to use the provided setters, but those are currently difficult to
use with Parameter optimization using ForwardDiff and Zygote.
In the meantime, MTKHelpers does not support non-scalar states or parameters.

See subpages of the [package help](https://bgctw.github.io/MTKHelpers.jl/dev)
