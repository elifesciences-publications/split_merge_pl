all:
	(cd others/BarcodeSplitter; make all)
	(cd others/read_trimmer; make all)
	(cd others/dual_indexed_bc; make all)
#	(cd shell_scripts; make all)
clean:
	(cd others/BarcodeSplitter; make clean)
	(cd others/read_trimmer; make clean)
	(cd others/dual_indexed_bc; make clean)
#	(cd shell_scripts; make clean)
