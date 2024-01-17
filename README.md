# PATSMA

## Running Example

Installing the PATSMA library

```sh
$ mkdir build && cd build
$ cmake ../ -DCMAKE_INSTALL_PREFIX=../example/patsma/ -DOPENMP=ON
$ make -j install
```

Installing the example

```sh
$ cd ../example
$ mkdir build && cd build
$ cmake ../ -DCMAKE_INSTALL_PREFIX=./
$ make -j install
```

Running the example

```sh
$ ./bin/insideat <dimension>
$ ./bin/outsideat <dimension>
$ ./bin/noat
```

## Citing PATSMA

```bibtex
@misc{patsma2024,
  title={PATSMA: Parameter Auto-tuning for Shared Memory Algorithms},
  author={Joao B. Fernandes and Felipe H. S. da Silva and Samuel Xavier-de-Souza and Italo A. S. Assis},
  year={2024},
  eprint={2401.07861},
  archivePrefix={arXiv},
  primaryClass={cs.DC}
}
```

```bibtex
@article{at-rtm-2020,
  author={Assis, Ítalo A. S. and Fernandes, João B. and Barros, Tiago and Xavier-De-Souza, Samuel},
  journal={IEEE Access},
  title={Auto-Tuning of Dynamic Scheduling Applied to 3D Reverse Time Migration on Multicore Systems},
  year={2020},
  volume={8},
  number={},
  pages={145115-145127},
  doi={10.1109/ACCESS.2020.3015045}
}
```