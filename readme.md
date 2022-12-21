# TEGBP
Code for CVPR submmision 4712.

- [Overleaf](https://www.overleaf.com/project/634f900abe6816625b0e860e)
- [PDF](https://github.com/DensoITLab/tegbp/blob/master/material/CVPR2023_4712.pdf)


<img src="https://github.com/DensoITLab/tegbp/blob/master/material/flo_1_00009.png" alt="normal flow"  title="normal flow">
<img src="https://github.com/DensoITLab/tegbp/blob/master/material/flo_0_00009.png" alt="full flow"  title="full flow">


https://user-images.githubusercontent.com/270088/208857729-874e8b45-a718-4e47-93a3-04094538ba18.mov


## Setup dependency 
```
sudo apt install libeigen3-dev
pip install Command
```

## Prepare result folder
Local SSD is recomended.
```
ln -s /home/data2/cashe/tegbp_result/ result
sudo chmod -R 777 result
mkdir result/bricks_1slide
mkdir result/dummy
mkdir result/bricks
mkdir result/indoor_flying2
```

## Compile 
```
make
```

## Run
```
./process #thread #events
./process =1 "indoor_flying2" -1 -1
./process -1 "bricks" 15000 // debug
```

## Visualization
Run `vis_result.ipynb`.
It'll execute the followng operations
- run the core program (cpp),
- show the optical flow,
- and saved the resutls as png file

If you can not open notebook, first convert it to the python script and then execute as a python script.
```
jupyter nbconvert  --to script vis_result.py  vis_result.ipynb 
python vis_result.py 
```

## Misc
See  `matlab` directory for the reference matlab implementation (slow).

## Credit
Coder: Y.Sekikawa, J.Nagata

## Reference
- [Compile option](https://ac-as.net/gcc-optimization-option/)
- [OMP basics](https://ichigopack.net/openmp/omp_base_3.html)
- [OMP hello](https://curc.readthedocs.io/en/latest/programming/OpenMP-C.html#id2)
- [OMP tutorial](https://msu-cmse-courses.github.io/cmse401-S21-student/assignments/0225-HW2_OMP.html)
- [OMP youtube](https://www.youtube.com/playlist?list=PLLX-Q6B8xqZ8n8bwjGdzBJ25X2utwnoEG)
