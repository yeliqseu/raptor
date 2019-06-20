Introduction
============
This is a prototype Raptor code implemention. It was intented for performance comparison with sparse network codes, see: <https://github.com/yeliqseu/sparsenc>. The implementation did not follow RFCs (5053, 6330) strictly. Instead, only some key features of R10 were implemented, including LDPC precoding, inactivation decoding, as well as the hardware acceleration of Galois field operations using SIMD (Intel SSSE3 and AVX2).

Usage
============
See test.c for a simulation example of using the code to send some bytes over a multi-hop lossylink. Compile it by

```shell
$ make raptorTest
```
