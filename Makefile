JL = julia --project

default: init test

init:
	$(JL) -e 'using Pkg; Pkg.precompile(); Pkg.activate("docs"); Pkg.develop(path="."); Pkg.precompile(); Pkg.activate("examples"); Pkg.develop(path="."); Pkg.precompile()'

update:
	$(JL) -e 'using Pkg; Pkg.update(); Pkg.precompile(); Pkg.activate("docs"); Pkg.update(); Pkg.precompile(); Pkg.activate("examples"); Pkg.update(); Pkg.precompile()'

test:
	$(JL) -e 'using Pkg; Pkg.test()'

coverage:
	$(JL) -e 'using Pkg; Pkg.test(; coverage=true)'

serve:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); using LiveServer; servedocs(;skip_dirs=["docs/src/assets", "docs/src/generated"], literate_dir="examples")'

showme-hypernd:  # QEC does not work with KaHyPar
	for case in inference quantumcircuit nqueens independentset; do \
		echo "Running $${case}"; \
		$(JL) -e "rootdir=\"examples/$${case}\"; using Pkg; Pkg.activate(rootdir); Pkg.develop(path=\".\"); Pkg.instantiate(); include(joinpath(rootdir, \"main.jl\")); using KaHyPar; main(HyperND())"; \
	done

showme-treesa:  # QEC does not work with KaHyPar
	for case in inference quantumcircuit nqueens independentset qec; do \
		echo "Running $${case}"; \
		$(JL) -e "rootdir=\"examples/$${case}\"; using Pkg; Pkg.activate(rootdir); Pkg.develop(path=\".\"); Pkg.instantiate(); include(joinpath(rootdir, \"main.jl\")); main(TreeSA())"; \
	done

update-examples:  # QEC does not work with KaHyPar
	for case in inference quantumcircuit nqueens independentset qec; do \
		echo "Running $${case}"; \
		$(JL) -e "rootdir=\"examples/$${case}\"; using Pkg; Pkg.activate(rootdir); Pkg.update()"; \
	done

fig:
	for entry in "docs/src/assets/"*.typ; do \
		echo compiling $$entry to $${entry%.typ}.pdf; \
		typst compile $$entry $${entry%.typ}.pdf; \
		pdf2svg $${entry%.typ}.pdf $${entry%.typ}.svg; \
	done

clean:
	rm -rf docs/build
	find . -name "*.cov" -type f -print0 | xargs -0 /bin/rm -f

.PHONY: init test coverage serve clean update