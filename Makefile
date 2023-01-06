.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<

## clean_data: remove all data and parameter files
.PHONY : clean_data
make clean_data:
	rm data/*.npz
	rm data/*.json
