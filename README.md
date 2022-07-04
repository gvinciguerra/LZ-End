# LZ-End parsing

This repository contains an implementation of two LZ-End parsing algorithms described in the papers:

- Sebastian Kreft and Gonzalo Navarro. [On compressing and indexing repetitive sequences](https://doi.org/10.1016/j.tcs.2012.02.006). Theor. Comput. Sci. (2013).
- Dominik Kempa and Dmitry Kosolobov. [LZ-End Parsing in Linear Time](https://doi.org/10.4230/LIPIcs.ESA.2017.53). ESA 2017. [Only the algorithm of §2]

The implementation supports both **byte** and **integer** alphabets.
It does not define a compressed file format, rather it outputs a list of LZ-End phrases, each represented as a `struct` with the phrase length, the id of a previous phrase where the phrase end, and the trailing symbol.

Check out also [this repository](https://github.com/dominikkempa/lz-end-toolkit) by Kempa and Kosolobov.

## Usage

This is a header-only library. To compile the [example](example.cpp), use the following commands:

```sh
git clone https://github.com/gvinciguerra/LZ-End.git
cd LZ-End
cmake . -DCMAKE_BUILD_TYPE=Release
make -j8
```

## License

This project is licensed under the terms of the Apache License 2.0.