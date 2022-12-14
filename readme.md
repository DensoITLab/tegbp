# TEGBP
Code for CVPR submmision 4712.

[Overleaf](https://www.overleaf.com/project/634f900abe6816625b0e860e)
[PDF](https://github.com/DensoITLab/tegbp/blob/master/material/CVPR2023_4712.pdf)

## Setup dependency 
```
sudo apt install libeigen3-dev
pip install Command
```


## Compile 
```
make
```

## Run
```
./process #thread #events
./process 64 100000
```

## Visualization
Run `vis_result.ipynb`.
It'll show the optical flow and the result are saved as png file.


## Misc
See  `matlab` directory for the reference matlab implementation (slow).

## Citation
Coder: Y.Sekikawa, J.Nagata

## Reference
- [Compile option](https://ac-as.net/gcc-optimization-option/)
- [OMP basics](https://ichigopack.net/openmp/omp_base_3.html)
- [OMP tutorial](https://msu-cmse-courses.github.io/cmse401-S21-student/assignments/0225-HW2_OMP.html)
