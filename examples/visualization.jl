using OMEinsumContractionOrders, LuxorGraphPlot

eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], ['a'])

viz_eins(eincode, filename = "eins.png")

nested_eins = optimize_code(eincode, uniformsize(eincode, 2), GreedyMethod())
viz_contraction(nested_eins)