using TensorQEC, Random, OMEinsumContractionOrders, CliqueTrees, Metis

d = 9
tanner = CSSTannerGraph(SurfaceCode(d, d))

optimizer = Treewidth(; alg=ND(MF(), METISND(; ufactor=150); limit=200, level=6));
ct = compile(TNMAP(; optimizer),tanner)

# NOTE: TreeSA gives sc = 25
contraction_complexity(ct.net)