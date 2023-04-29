

BASE_URL="https://storage.googleapis.com/jlumpe-gambit/public/publication/genomes/set3/fasta"
GENOMES="16AC1611138-CAP,17AC0001410A,17AC0006310,17AC0006313-1,19AC0011210"

.PHONY: all


all: results.csv results-transposed.csv example-genomes.zip


# Download genome files
genomes:
	mkdir -p "genomes"
	curl --silent -C - --remote-name-all "${BASE_URL}/{${GENOMES}}.fasta.gz" -o "genomes/#1.fasta.gz"


# Gambit query output
results.csv: genomes
ifndef GAMBIT_DB_PATH
	$(error GAMBIT_DB_PATH environment variable must be set)
endif
	gambit query -o $@ $</*


# Transposed table for easier reading
results-transposed.csv: results.csv
	python scripts/transpose-results.py $@ $<


# Put genomes in zip file for download
example-genomes.zip: genomes
	zip -0 $@ genomes/*
