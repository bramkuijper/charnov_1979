# Individual-based simulation of Charnov (1979)
Evolution of sex allocation in hermaphrodites with ova fitness wf=r and pollen fitness `wm=r^n`. Written in C++.

See Charnov (1979) PNAS [https://doi.org/10.1073/pnas.76.5.2480](https://doi.org/10.1073/pnas.76.5.2480)

## How to run?
On any bash shell:
```cd src/ibm
make
./sexalloc.exe
```
Resulting data will be retained in the file `data.csv` which can be plotted using the corresponding `plot_sexalloc.r` R script.

## Examples
Charnov's candidate ES predicts `r* = n/(n+1)`. Hence, for `n=0.3` we have `r*=0.23` which is pretty close to what the simulation shows in the figure below:

![n03](https://github.com/bramkuijper/charnov_1979/blob/main/img/rplot_n03.png?raw=true)

Similarly, for `n=0.8` we would predict `r*=0.45` again close to what we find here:

![n08](https://github.com/bramkuijper/charnov_1979/blob/main/img/rplot_n08.png?raw=true)

For values of `n>1.0` we find evolutionary branching into males and females as expected.
