from mersenne_twister import MersenneTwister
from z3 import *


class MersenneCracker:
    """[summary]
    https://nayak.io/posts/mersenne_twister/

    观察到的随机数序列需要大于或者等于n

    有两种使用方式：
    1. 已知连续的n个随机数且没有发生twist -> 使用crack_state (注意代码中是取前n个随机数)
    2. 已知连续的n个随机数但不确定twist的情况 -> 使用_brute_force_closest_twist_index或者predict_next_state
    """    
    def __init__(self, variant = "mt19937", parameters = []):
        if variant == "mt19937":  # 32-bit version of original twister
            self.w = 32
            self.n = 624
            self.m = 397
            self.r = 31
            self.a = 0x9908b0df
            self.u = 11
            self.d = 0xffffffff
            self.s = 7
            self.b = 0x9d2c5680
            self.t = 15
            self.c = 0xefc60000
            self.l = 18
            self.f = 1812433253
        elif variant == "mt19937_64":  # 64-bit version of original twister
            self.w = 64
            self.n = 312
            self.m = 156
            self.r = 31
            self.a = 0xb5026f5aa96619e9
            self.u = 29
            self.d = 0x5555555555555555
            self.s = 17
            self.b = 0x71d67fffeda60000
            self.t = 37
            self.c = 0xfff7eee000000000
            self.l = 43
            self.f = 6364136223846793005
        elif variant == "mt11213b":  # mt11213b twister version from boost
            self.w = 32
            self.n = 351
            self.m = 175
            self.r = 19
            self.a = 0xccab8ee7
            self.u = 11
            self.d = 0xffffffff
            self.s = 7
            self.b = 0x31b6ab00
            self.t = 15
            self.c = 0xffe50000
            self.l = 17
            self.f = 1812433253
        elif variant == "custom":  # only use if you know what you're doing
            assert len(parameters) == 13
            self.w = parameters[0]
            self.n = parameters[1]
            self.m = parameters[2]
            self.r = parameters[3]
            self.a = parameters[4]
            self.u = parameters[5]
            self.d = parameters[6]
            self.s = parameters[7]
            self.b = parameters[8]
            self.t = parameters[9]
            self.c = parameters[10]
            self.l = parameters[11]
            self.f = parameters[12]
        else:
            raise NotImplementedError
        self.lower_mask = int("1"*(self.r), 2)
        self.upper_mask = int("1"*(self.w-self.r)+"0"*self.r, 2)


    def crack_state(self, outputs):
        self.original_state = [0] * self.n
        assert len(outputs) >= self.n
        assert all(isinstance(i, int) for i in outputs)
        outputs = outputs[:self.n]

        # reverses the temper operations
        self.original_state = map(self._untemper, outputs)
        return list(self.original_state)

    def untemper_right(self, n, shift):
        i = 0
        while i * shift < self.w:
            new_mask = n & (((((1 << self.w) - 1) << (self.w - shift)) & ((1 << self.w) - 1)) >> (shift * i))
            n ^= new_mask >> shift
            i += 1
        return n

    def untemper_right_mask(self, n, shift, mask):
        i = 0
        while i * shift < self.w:
            new_mask = n & (((((1 << self.w) - 1) << (self.w - shift)) & ((1 << self.w) - 1)) >> (shift * i))
            new_mask >>= shift
            n ^= new_mask & mask
            i += 1
        return n

    def untemper_left(self, n, shift, mask):
        i = 0
        while i * shift < self.w:
            new_mask = n & ((((1 << self.w) - 1) >> (self.w - shift)) << (shift * i))
            new_mask <<= shift
            n ^= new_mask & mask
            i += 1
        return n

    def untwist(self, state):  # returns the past state before a twist
        state = list(state)
        for i in range(self.n - 1, -1, -1):
            temp = state[i] ^ state[(i + self.m) % self.n]  # find leading bit
            if temp % 2:
                temp ^= self.a
            shifted = (temp << 1) & self.upper_mask
            temp = state[(i - 1) % self.n] ^ state[(i + self.m - 1) % self.n]  # ending bits of int
            if temp & self.upper_mask == self.upper_mask:  # check if leading bit is the same
                temp ^= self.a
                shifted |= 1
            state[i] = shifted ^ (temp << 1) & self.lower_mask
        return state


    def untemper(self, observed_values):
        """returns untempered values of observed
        """        
        return list(map(self._untemper, observed_values))

    def _untemper(self, single_observed_value):
        x = single_observed_value
        y = self.untemper_right(x, self.l)
        y = self.untemper_left(y, self.t, self.c)
        y = self.untemper_left(y, self.s, self.b)
        y = self.untemper_right_mask(y, self.u, self.d)
        return y

    def _untemper_with_z3(self, single_observed_value):
        """https://www.schutzwerk.com/en/43/posts/attacking_a_random_number_generator/
        x = self.state[self.index]
        y1 = x ^ ((x >> self.u) & self.d) 
        y2 = y1 ^ ((y1 << self.s) & self.b)
        y3 = y2 ^ ((y2 << self.t) & self.c)
        y ^= (y3 >> self.l)

        特别注意: 
        1. 需要使用逻辑右移(移走的位填充为0)
        2. 需要使用算术左移(尾部补0)
        我们就是取高位和低位, 所以无论如何都是填0

        https://docs.python.org/3/reference/expressions.html#shifting-operations
        Python only lets you do the arithmetic shift
        A right shift by n bits is defined as floor division by pow(2,n).
        A left shift by n bits is defined as multiplication with pow(2,n).

        对于正数, 逻辑和算术是一样的结果, 故前面的temper不影响
        但在这里我们指定了位数, 故实现的时候需要区别

        https://realpython.com/python-bitwise-operators/#arithmetic-vs-logical-shift
        """
        x = BitVec('x', self.w)
        y1 = BitVec('y1', self.w)
        y2 = BitVec('y2', self.w)
        y3 = BitVec('y3', self.w)
        y = BitVecVal(single_observed_value, self.w)
        s = Solver()
        equations = [
            y1 == x ^ (LShR(x, self.u) & self.d),
            y2 == y1 ^ ((y1 << self.s) & self.b),
            y3 == y2 ^ ((y2 << self.t) & self.c),
            y == y3 ^ LShR(y3, self.l)
        ]
        
        s.add(equations)
        assert s.check() == sat
        return s.model()[x].as_long()
    

    def _predict_next_state(self, x, closest_twist_index):
        assert len(x) >= self.n
        second_closest_twist_index = closest_twist_index - self.n
        k_plus_n = len(x) # next state
        k = k_plus_n -self.n
        k_plus_1 = (k_plus_n -self.n + 1 - closest_twist_index) % self.n + second_closest_twist_index
        k_plus_m = (k_plus_n -self.n + self.m - closest_twist_index) % self.n + second_closest_twist_index
        x__k = x[k]
        x__k_plus_1 = x[k_plus_1]
        x__k_plus_m = x[k_plus_m]

        # twist: x[k+n] only depends on x[k], x[k+m], x[k+l]
        temp = (x__k & self.upper_mask) + (x__k_plus_1 & self.lower_mask)
        shifted = temp >> 1
        if temp % 2: # if the lowest order bit of x==1 -> xor self.a
            shifted ^= self.a
        x__k_plus_n = x__k_plus_m ^ shifted
        return self._temper(x__k_plus_n)


    def _brute_force_closest_twist_index(self, observed_values, x):
        """[summary]

        x[k+n] only depends on x[k], x[k+m], x[k+l]

        Assume observed_values[:- (len(observed_values) - self.n)] as unknown, brute force twist index
        """
        assert len(observed_values) > self.n
        assert len(observed_values) == len(x)

        best_guesses = [] # brute force twist index
        best_guess_right_count = 0

        for tmp_twist_index in range(self.n, len(observed_values)):
            
            tmp_best_guess_right_count = 0
            for guess_index in range(0, len(observed_values) - self.n - 1):
                try:
                    y_hat = self._predict_next_state(x[:-1-guess_index], tmp_twist_index)
                    # print(tmp_twist_index,guess_index, y_hat, observed_values[-1-guess_index])
                    if y_hat == observed_values[-1-guess_index]:
                        tmp_best_guess_right_count += 1
                    else:
                        break
                except:
                    pass

            if tmp_best_guess_right_count > best_guess_right_count:
                best_guess_right_count = tmp_best_guess_right_count
                best_guesses.clear()
                best_guesses.append(tmp_twist_index)
            elif tmp_best_guess_right_count == best_guess_right_count:
                best_guesses.append(tmp_twist_index)
        
        return list(best_guesses)

    def _temper(self, x):
        """Tempered representation of internal states
        where x is the next value from the series, 
        """        
        y = x ^ ((x >> self.u) & self.d) 
        y ^= ((y << self.s) & self.b)
        y ^= ((y << self.t) & self.c)
        y ^= (y >> self.l)
        return self._fixed_int(y) # lowest w bits of y
 
    def _fixed_int(self, n):  # 
        """truncate int to lowest w bits
        ((1 << self.w) - 1) = 0b11111111111111111111111111111111 if w=32
        """        
        return ((1 << self.w) - 1) & n

    def predict_next_state(self, observed_values):
        predictions = set()
        x = self.untemper(observed_values) # states
        guessed_closest_twist_index_list = self._brute_force_closest_twist_index(observed_values, x)
        for guessed_closest_twist_index in guessed_closest_twist_index_list:
            predictions.add(self._predict_next_state(x, guessed_closest_twist_index))
        return list(predictions)

if __name__ == "__main__":
    # 第一种使用场景: 已知连续的n个随机数且没有发生twist -> 使用crack_state
    # 对于随机数个数超过n的情况, 默认是取前n个 (line 83)
    random_64 = MersenneTwister(variant = "mt19937_64")
    outputs_64 = [random_64.random_integer() for _ in range(1024)]
    cracker_64 = MersenneCracker(variant = "mt19937_64")

    # 假设我们只观察到了前n个随机数, 后面的随机数未观察到
    observed_values = outputs_64[:random_64.n]
    new_state_64 = cracker_64.crack_state(observed_values)
    random_64.setstate([new_state_64, 0])  # start at beginning of state
    assert outputs_64 == [random_64.random_integer() for _ in range(1024)]
    print("64-bit Case 1 Successfully Cracked")

    # ===================================================
    # 第二种使用场景: 已知连续的n个随机数但不确定twist的情况 -> 使用_brute_force_closest_twist_index或者predict_next_state
    observed_values = []
    with open("520.txt", "r") as f:
        for line in f.readlines():
            num = int(line)
            observed_values.append(num)
    cracker_64 = MersenneCracker(variant = "mt19937_64")
    guessed_random_numbers = cracker_64.predict_next_state(observed_values)
    assert 310574417605022995 in guessed_random_numbers
    print("64-bit Case 2 Successfully Cracked")

