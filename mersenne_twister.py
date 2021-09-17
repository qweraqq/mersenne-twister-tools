class MersenneTwister:
    """梅森旋转算法
    https://en.wikipedia.org/wiki/Mersenne_Twister

    w: word size (in number of bits)
    n: degree of recurrence
    m: middle word, an offset used in the recurrence relation defining the series x, 1 ≤ m < n
    r: separation point of one word, or the number of bits of the lower bitmask, 0 ≤ r ≤ w − 1
    a: coefficients of the rational normal form twist matrix
    b, c: TGFSR(R) tempering bitmasks
    s, t: TGFSR(R) tempering bit shifts
    u, d, l: additional Mersenne Twister tempering bit shifts/masks

    """    
    def __init__(self, mt_seed = 5489, variant = "mt19937", parameters = []):
        mt_seed = hash(mt_seed)  # hashes non-int types into ints, ints are unchanged
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
        self.lower_mask = (1 << self.r) - 1 # 0b1111111111111111111111111111111 if r=31
        self.upper_mask = 1 << self.r # 0b10000000000000000000000000000000 if r=31
        self.seed(mt_seed) # Initialization

    def random_integer(self):
        """Tempered representation of internal states
        where x is the next value from the series, 
        y a temporary intermediate value, 
        z the value returned from the algorithm, 
        with ≪, ≫ as the bitwise left and right shifts, 
        and & as the bitwise and. 
        
        The first and last transforms are added in order to improve lower-bit equidistribution. 
        From the property of TGFSR, s + t ≥ ⌊w/2⌋ − 1 is required to reach the upper bound of equidistribution for the upper bits.
        """        
        if self.index >= self.n:
            self.twist()
        x = self.state[self.index]
        y = x ^ ((x >> self.u) & self.d) 
        y ^= ((y << self.s) & self.b)
        y ^= ((y << self.t) & self.c)
        y ^= (y >> self.l)
        self.index += 1
        return self.fixed_int(y) # lowest w bits of y

    def twist(self):
        """Generate the next n values from the series x_i
        x[k+n] only depends on x[k], x[k+m], x[k+l]

        也就是说下个x_k的状态仅依赖于现在的x[k] x[k+m] x[x+l]
        """        
        for i in range(self.n):
            temp = (self.state[i] & self.upper_mask) + (self.state[(i + 1) % self.n] & self.lower_mask)
            shifted = temp >> 1
            if temp % 2:
                shifted ^= self.a
            self.state[i] = self.state[(i + self.m) % self.n] ^ shifted
        self.index = 0

    def seed(self, seed):
        """Initialization
        The state needed for a Mersenne Twister implementation is an array of n values of w bits each. 
        To initialize the array, a w-bit seed value is used to supply x[0] through x[n−1] by setting x[0] to the seed value and thereafter setting

        x[i] = some_function( x[i-1] ) for i to n-1
        
        Args:
            seed ([type]): [description]
        """        
        self.state = [0] * self.n
        self.state[0] = seed
        self.index = self.n
        for i in range(1, self.n):
            self.state[i] = self.fixed_int(self.f * (self.state[i - 1] ^ (self.state[i - 1] >> (self.w - 2))) + i)

    def getstate(self):  # [state, index] format
        return [list(self.state), self.index]

    def setstate(self, new_state):  # [state, index] format
        assert len(new_state[0]) == self.n
        assert all(isinstance(i, int) for i in new_state[0])
        self.state = new_state[0]
        self.index = new_state[1]

    def fixed_int(self, n):  # 
        """truncate int to lowest w bits
        ((1 << self.w) - 1) = 0b11111111111111111111111111111111 if w=32
        """        
        return ((1 << self.w) - 1) & n

if __name__ == "__main__":
    random_32 = MersenneTwister(mt_seed=123)
    print(f"Test 32-bit {[random_32.random_integer() for i in range(32)]}")

    random_64 = MersenneTwister(variant = "mt19937_64")  # 64-bit numbers
    print(f"Test 64-bit: {random_64.random_integer()}")
    random_alt = MersenneTwister(variant = "mt11213b")  # alternate twister
    print(f"Test Alt 32-bit: {random_alt.random_integer()}")
    saved_state = random_32.getstate()
    print(f"Original State: {random_32.random_integer()}")  # advances the state
    random_32.setstate(saved_state) # restores the state
    print(f"Restored State: {random_32.random_integer()}")  # verifies the restored state
