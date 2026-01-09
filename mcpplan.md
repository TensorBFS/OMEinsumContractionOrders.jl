# MCP Helper

1. Create a package named `MCPHelper.jl`.
2. This package contains a macro `@mcp` that takes a function as input, and
   - make sure all input arguments are simple variables and can be directly translated to JSON.
   - make sure the output is asserted to be a named tuple.
   - returns an MCP tool.

Ref: the function expression analysis and meta programming follow the style of the `@register` macro of the `Jugsaw` package.