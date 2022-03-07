# Differential expression

## Folder organization

..* src (This is where the functions reside)
..* data (This is where data reside)
..* report (This is where the code and pipelines reside)

### Prerequisites. R version 3.4.1

The following packages should to be installed:
⋅⋅* "gplots"
⋅⋅* "tidyverse"
⋅⋅* "broom"
⋅⋅* "here"
⋅⋅* "data.table"
⋅⋅* "RUVSeq"
⋅⋅* "edgeR"
⋅⋅* "statmod"
⋅⋅* "RColorBrewer"
⋅⋅* "drake"

### Data
Raw RNA count data is found in /data/counts52KPS.txt

### Example
```r
install.packages("gplots")
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
```

### Building the project

The following should be followed to obtain the data inside of R:
1. Open the report folder
2. Start R(version 3.4.1)
3. Run source(make.R)
4. The variables are now loaded.

### Example
```r
dat <- readd(dat)
DE-results <- readd(FDRres)
```

You can now work with dat (The raw count data)

## Authors

* **Thomas Wolfe** - *Initial work* - [twolfe](https://github.com/twolfe)

See also the project website of [contributors](http://www.botanik.univie.ac.at/systematik/projects/dactylorhiza/people.html) who participated in this project.

## License

Users are free to use and modify code.
