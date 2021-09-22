# 伪随机数生成器
## 真随机与伪随机
- 自然 (non-deterministic): 如从电子噪声采样
- machines (deterministic)
## 早期的PRNG (Enigma )
- by JOHN VON NEUMAN
- Middle squares method
```
seed = 43824
square: 43824 ^ 2 = 1920542976
output the middle: 92054297

next seed = 92054297
...
```

------
# PRNG vuls
## Choose a "secure" seed
- time
- pids
- stack memory
- hardcoded random seed
- comninations of any of the above
## Attacking
### Brute force
- 32 bit seeds; offline attack; timestamp based seeds
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

------
# Mersenne twister
## About Mersenne twister
- PRNG period size: `2^19937 -1` 
- [wikipedia](https://en.wikipedia.org/wiki/Mersenne_Twister)
```
Is not cryptographically secure, unless the CryptMT variant (discussed below) is used. The reason is that observing a sufficient number of iterations (624 in the case of MT19937, since this is the size of the state vector from which future iterations are produced) allows one to predict all future iterations.
```
## Attacks
MersenneCracker有两种使用方式：
1. 已知连续的n个随机数且没有发生twist -> 使用`crack_state` (注意代码中是取前n个随机数)
2. 已知连续的n个随机数但不确定twist的情况 -> 使用`_brute_force_closest_twist_index`或者`predict_next_state`

------
# Credits
- [https://en.wikipedia.org/wiki/Mersenne_Twister](https://en.wikipedia.org/wiki/Mersenne_Twister)
- [Bsides LV 2014 - Untwisting The Mersenne Twister: How I killed the PRNG - 05Aug2014](https://www.youtube.com/watch?v=f841Y7d3oDo&list=PLcO_vga2cLLu0PFSrJ8nsO92YBdC2lqwm&index=1)
- [https://nayak.io/posts/mersenne_twister/](https://nayak.io/posts/mersenne_twister/)
- [https://www.schutzwerk.com/en/43/posts/attacking_a_random_number_generator/](https://www.schutzwerk.com/en/43/posts/attacking_a_random_number_generator/)
- [https://realpython.com/python-bitwise-operators/#arithmetic-vs-logical-shift](https://realpython.com/python-bitwise-operators/#arithmetic-vs-logical-shift)
- [https://github.com/slightlyskepticalpotat/mersenne-twister-tools](https://github.com/slightlyskepticalpotat/mersenne-twister-tools)
