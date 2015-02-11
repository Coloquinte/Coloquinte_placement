# Coloquinte: an analytical VLSI placer

Coloquinte is a framework for placing VLSI circuit. It targets standard-cell digital circuits, optionally with pre-placed macroblocks.

It is now a complete tool, handling all steps of a placement flow. It features exact and approximate wirelength models, and should be updated to do timing- and routing-driven optimization.

It is not standalone: no parsers; currently, it interfaces with the Coriolis CAD toolchain, but you you can follow the example main and interface it with your own tool.
It is still a work in progress, and the optimization process could be significantly improved in the coming months. However, all tools are already competitive with other academic placers - and it is the only open source placer using analytical placement techniques.

## Principle

Coloquinte is based on analytical placement: it performs continuous optimization of the objective function with a penalty accounting for overlap. This penalty is currently calculated from a legalized or partially legalized placement.
It alternates between continuous optimization and legalization to reach a good solution, then perform local optimizations.

## Tools

### Global placement

Currently, the global placement is based on a common method: a quadratic local wirelength model that is optimized using conjugate gradient.
It is extremely flexible, with the exact Steiner model as well as the classical star and bounding-box models, and takes advantage of parallelism.

### Legalization

Legalization is performed in two steps: rough legalization, accounting for placement density, and exact legalization accounting for cell overlaps.
Both should be beyond industrial level, although exact legalization only handles standard cells yet.

The rough legalizer models a 2D transportation problem that is optimized using various locally optimal algorithms. They may be run in parallel.
The exact legalization places one cell at a time but is able to move previous cells to obtain better solution: it can handle extremely dense placements (>99%) with fixed macroblocks.

### Detailed placement

The detailed placement optimizations look for topology (i.e. cell ordering) modifications, and optimize the positions for a given topology.
They use specialized algorithms for position optimizations - most of them are not published yet -, and generally brute-force for topology modifications.

## Runtime

For ~300,000 cells and a single core at 2 GHz, the runtimes are approximately (depending on solution quality):
  * 10-30s per global placement iteration
  * 10s for a rough legalization
  * 3s to obtain a fully legal placement
  * 400s to perform detailed placement

A typical placement run will usually feature 50 to 200 global placement iterations. Note that global placement and rough legalization can use multiple cores.








