# cp2kmd2g96

A CP2K AIMD output to g96 convertor

## Build

```sh
$ make
```

or

```sh
$ FC=ifort make
```

## Usage

```sh
$ ./cp2kmd2g96.x
Usage:
  $ cp2kmd2g96.x CELL POS [VEL] [G96]
```

```sh
$ ./cp2kmd2g96.x md-1.cell md-pos-1.xyz traj.g96
$ ./cp2kmd2g96.x md-1.cell md-pos-1.xyz md-vel-1.xyz traj.g96
```

## Example

See `example` folder.

## License

BSD 2-Clause License
