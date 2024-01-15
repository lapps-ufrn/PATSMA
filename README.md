# PATSMA

Library documentation in https://lappsufrn.gitlab.io/auto-tuning/.

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