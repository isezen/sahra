# sahra
Saharan Dust Transport Research

##### How to install sahra to your computer?

Run the command below in your terminal or use GitHub Desktop for OSX users.

```
git clone https://github.com/isezen/sahra.git
```

##### Prerequisities

Install packages below before use the sahra.
```R
devtools::install_github("isezen/rwind")
devtools::install_github("isezen/rsezen")
```

##### setup.sh

Run `setup.sh` after you downloaded/cloned sahra to your computer. `setup.sh` downloads required netcdf files to your `data` directory and converts them to `rdata` format for easy processing.

```bash
./setup.sh
```
