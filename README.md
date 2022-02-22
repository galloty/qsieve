# qsieve
A Quad Sieve for the search for twin and Sophie Germain primes

## About

Search for quadruples (*k* 2<sup>*n* - 1</sup> - 1, *k* 2<sup>*n*</sup> - 1, *k* 2<sup>*n*+1</sup> - 1, *k* 2<sup>*n*</sup> + 1) such that neither of them has a small factor.

We must have *k* = 15 (mod 30). If the sieved range is *k*<sub>*min*</sub> = 3, *k*<sub>*max*</sub> = 2<sup>32</sup>, *n*<sub>*min*</sub> = 3321925 and *n*<sub>*max*</sub> = *n*<sub>*min*</sub> + 63 then the bitmap size is 1 gigabyte.

*qsieve* is multithreaded:
- Thread 1 is a pseudo prime number generator based on a segmented sieve of Eratosthenes.
- Thread 2 computes 1 / 15 · 2<sup>*n*_*min* - 1</sup> mod *p* for each *p*.
- Thread 3 fills the bitmap for all *k* and *n* in the range such that *k* 2<sup>*n*</sup> &plusmn; 1 = 0 (mod *p*).

Δ*n* is set to 64. This enables the bitmap to be a fast 64-bit array.

The expected number of remaining candidates is 0.41252 · Δ*k* · Δ*n* / (log *p*<sub>*max*</sub>)<sup>4</sup>.<br>
For 3 &le; *k* < 2<sup>32</sup>, 3321925 &le; *n* &le; 3321988 and *p*<sub>*max*</sub> = 10<sup>12</sup>, 194534 candidates are expected and their actual number is 195115.

## Build

This version was compiled with gcc 11 and C++17 standard.
