# 伪随机数生成器

## 真随机与伪随机
- 自然 (non-deterministic): 如从电子噪声采样
- machines (deterministic)

## 早期的PRNG
- Enigma 
- by JOHN VON NEUMAN
- Middle squares method
```
seed = 43824
square: 43824 ^ 2 = 1920542976
output the middle: 92054297

next seed = 92054297
...
```
# PRNG vuls 
## Choose a "secure" seed
- time
- pids
- stack memory
- hardcoded random seed
- comninations of any of the above
## Attacking
### Brute force
- 32 bit seeds; o ffline attack; timestamp based seeds
- Algorithm
1. Given set of observed numbers (may not be fully observed)
2. Calculate random numbers and compare
    * If matches, look for the next observed
    * Continue until match is found
3. Allow for missing values

### State inference
- Example:  glibc rand() produces 31 bit numbers (with internal state has 32bit numbers)
If we observe 32 sequential numbers -> we get internal states (missing the LSBs); LSB reconstruction: guess and check
```
r[i] = r[i-3] + r[i-31]
with the LSB chopped off
```


# Mersenne twister
- PRNG period size: `2^19937 -1` 
- [wikipedia](https://en.wikipedia.org/wiki/Mersenne_Twister)
```
Is not cryptographically secure, unless the CryptMT variant (discussed below) is used. The reason is that observing a sufficient number of iterations (624 in the case of MT19937, since this is the size of the state vector from which future iterations are produced) allows one to predict all future iterations.
```


# mersenne-twister-tools
A collection of various programs implementing and/or related to the Mersenne Twister PRNG. This is intended to be a bare-minimum Python 3 implementation that users can build off of.

## Features
- Seeding with any hashable type*
- Saving and restoring the state
- Pseudo-random ints and floats
- Can be used to subclass [Random](https://docs.python.org/3/library/random.html)
- Derive past generated numbers

- Original 32-bit Mersenne Twister
- Original 64-bit Mersenne Twister
- Boost's mt11213b Mersenne Twister
- Custom Parameters Mersenne Twister

- Original 32-bit Mersenne Twister Cracker
- Original 64-bit Mersenne Twister Cracker
- Boost's mt11213b Mersenne Twister Cracker
- Custom Parameters Mersenne Twister Cracker

\*Set the environment variable `PYTHONHASHSEED = 0` to fix hashes for non-ints.

## Contributing
Please report any bugs that you encounter, and open a pull request if you have something to add to the program. Thank you!

## References
When I was creating this program, I found several resources that are useful for anyone seeking to understand the Mersenne Twister. They are listed here for reference:
- https://en.wikipedia.org/wiki/Mersenne_Twister
- http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/emt.html
- https://sci-hub.se/10.1145/272991.272995
- http://www.quadibloc.com/crypto/co4814.htm
- https://jazzy.id.au/2010/09/22/cracking_random_number_generators_part_3.html
