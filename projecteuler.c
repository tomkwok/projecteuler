#include <stdio.h> // printf(), FILE, fopen("file", "mode"), fgets(ln, lnlen, fp), fclose(fp), sscanf(src, "format", dest)
#include <stdlib.h> // malloc(size), calloc(size, sizeof(...))
#include <string.h> // strcat(dest, src), memcpy(dest, src, size)
#include <math.h> // sqrt(), pow(), log(), log10(), floor(), ceil()
#include <stdint.h>
#include <stdbool.h> // bool, true, false // GCC: 4-bit / 1-byte

typedef uint8_t uchar; // GCC: 4-bit / 1-byte (Note: abbreviation is made up by myself)
//typedef uint16_t ushort; // GCC: 16-bit / 2-byte (short)
typedef uint32_t uint; // GCC: 32-bit / 4-byte (int)
typedef uint64_t ulong; // GCC: 64-bit / 8-byte (long int)

/* Note for arrays: size + 1 for the null terminator */

#define string(len) ( calloc(len + 1, sizeof(char)) )
#define uintArray(size) ( malloc((size + 1) * sizeof(uint)) )

#define uintBit (8 * sizeof(uint)) // 32
#define uintArrayBitCleared(size) ( calloc(ceil(size / uintBit) + 1, sizeof(uint)) )
#define setBit(a, k) ( a[k / uintBit] |= 1 << k % uintBit ) // k >> 5, k & 31
#define clearBit(a, k) ( a[k / uintBit] &= ~(1 << k % uintBit) )
#define getBit(a, k) ( a[k / uintBit] & 1 << k % uintBit )

#define ucharArray(size) ( malloc((size + 1) * sizeof(uchar)) )

#define ucharArraySquare(side) ( malloc((side * side + 1) * sizeof(uchar)) )
#define getCell(a, side, row, col) ( a[row * side + col] )
#define ucharArrayTriangle(row) ( malloc((sum(row) + 1) * sizeof(uchar)) )

#define charToInt(c) ( c - '0' )
//#define max(a, b) ( a > b ? a : b )

// 1. sum of multiples of a or b inclusive

uint q001 (const uint a, const uint b, const uint n) {
	uint i, sum = 0;
	for (i = 0; i < n; i++)
		if (!(i % a) || !(i % b)) // divisible by a or b inclusive
			sum += i;
	return sum;
}

uint sum (const uint n) {
	return n * (n + 1) / 2; // Note: cannot be re-ordered, or .5 will be rounded off (floor'd) to integers
}

uint sum_div_by (const uint n, const uint a) {
	return a * sum(n / a);
}

uint q001_fast (const uint a, const uint b, uint n) {
	n--; // under n, as in the question, i.e. max. = n - 1
	return sum_div_by(n, a) + sum_div_by(n, b) - sum_div_by(n, a * b); // or multiples of (a * b) will be double counted
}

// 2. sum of even Fibonacci numbers

uint q002 (const uint n) {
	uint sum = 0, fo = 1, f = 1; // F old, F
	while (f < n) {
		if(!(f % 2))
			sum += f;
		uint fn; // F new
		fn = f + fo;
		fo = f;
		f = fn;
	}
	return sum;
}

uint q002_alt (const uint n) { // avoid modulo
	uint sum = 0, a = 1, b = 1, c = a + b; // compilation time optimization
	while (c < n) {
		sum += c; // (3n)th number in Fibonacci sequence
		a = c + b;
		b = a + c;
		c = b + a;
	}
	return sum;
}

uint q002_fast (const uint n) { // even Fibonacci numbers form another sequence
	uint sum = 0, fo = 0, f = 2; // even: **0**, 1, 1, [2], 3, 5, [8], ... // even: **0** *(imaginary Fibonacci number)*, 1, 1,
	while (f < n) {
		sum += f;
		uint fn; // F new
		fn = 4 * f + fo; // F = F(n-1) + F(n-2) = 4 * F(n-3) + F(n-6)
		fo = f;
		f = fn;
	}
	return sum;
}

// 3. largest prime factor

ulong q003 (ulong n) { // prime factorization
	// this function effectively is primality check (q003r(n) = n indicates n is prime)
	uint p = 2; // the first prime
	bool q = 1; // 6k +/- 1 flag, set to true/+1 first [alternative: int q = 1;]
	while (p <= sqrt(n)) {
		if (!(n % p)) {
			return q003(n / p); // recursion, optimized by GCC?
			//n /= p;
			//p = 2; // reset
			//q = 1; // reset
			continue;
		}
		switch (p) { // find the next probable prime
			case 2:
				p = 3;
				break;
			case 3:
				p = 5;
				break;
			/* prime after 2, 3, [5] must be in the form of (6k +- 1) since
			 *  6k is divisible by 2, 3, 6;
			 *  6k +- 2 is divisible by 2;
			 *  6k + 3 is divisible by 3.
			 * even if p is not prime, it is still coprime with n (recursion: n / p)
			 */
			default:
				p += (q ? 2 : 4); // alternative: p += 3 - q;
				q = !q; // alternative: q = -q;
		}
	}
	return n;
}

// 4. largest (even-digit) palindrome product of two n-digit numbers

/* Note: single-digit number is a palindrome by definition */

ulong reverse_ulong (ulong n) {
	ulong r = 0;
	while (n) {
		r *= 10; // r might overflow when reversed number is larger than the original
		r += n % 10;
		n /= 10;
	}
	return r;
}

bool is_palindrome (ulong n) { // avoid overflowing, but slow
	ulong div = 1;
	while (n / div >= 10)
		div *= 10;
	while (n) {
		int l = n / div;
		int r = n % 10;
		if (l != r)
			return 0;
		n = (n % div) / 10; // remove the leftmost and the rightmost digit
		div /= 100;
	}
	return 1;
}

/* A number is divisible by 11 if
 * (sum of odd digits - sum of even digits) is divisible by 11,
 * which is the case of any even-digit palindrome.
 * Assume the palindrome product start with 9, according to a multiplication table,
 * try only {1,9}, {3,3}, {7,7} to get 9 as the last digit of the product.
 */

ulong q004 (const uint d) { // heavily optimized, to be removed
	// uint c = 0; // iteration counter
	ulong i, j, p, pmax = 0, maxj;
	uint stepj;
	const ulong min = 9 * pow(10, d - 1); // for both i-loop and j-loop
	const ulong maxi = pow(10, d) - 1;
	
	for (i = maxi; i > min; i -= 2) {
		if (i % 11) { // j is divisible by 11 (5/6 probability in this loop)
			maxj = i - i % 11; // alternative: maxj = i / 11 * 11;
			stepj = 11;
		} else { // i is divisible by 11 (1/6 probability in this loop)
			maxj = i;
			stepj = 2;
		}
		for (j = maxj; j > min; j -= stepj) {
			// c++;
			p = i * j;
			if (p < pmax)
				break;
			if (p == reverse_ulong(p))
				pmax = p;
		}
	}
	// printf("q004 c=%d\n", c);
	return pmax;
}

ulong q004_fast (const uint d) {
	uint c = 0; // iteration counter
	ulong i, j, p, pmax = 0;
	const ulong pmin = pow(10, d * 2) - pow(10, d * 2 - d / 2); // 990,000
	const ulong min = pow(10, d) - pow(10, d - d / 2); // 900
	const ulong max = pow(10, d) - 1; // 999
	uint i10, j10;
	
	for (i = max - max % 11; i > min; i -= 11) {
		c++;
		if (!(i % 2)) // Skip 2, 4, 6, 8, 0 in the last digit, i.e. even numbers
			continue;
		i10 = i % 10;
		switch (i10) {
			case 5: // Skip 5 in the last digit, i.e. multiples of 5
				continue;
			case 9:
				j10 = 1;
				break;
			case 1:
				j10 = 9;
				break;
			default: // 3, 7
				j10 = i10;
		}
		for (j = max - (9 - j10); /* j > min */ ; j -= 10) {
			c++;
			p = i * j;
			//printf("q004_fast p=%ld, i=%ld, j=%ld\n", p, i, j);
			if (p < pmin || p < pmax)
				break;
			if (p == reverse_ulong(p))
				pmax = p;
		}
	}
	//printf("q004_fast d=%d c=%d\n", d, c);
	return pmax;
}

// 5. LCM of 1 to n

uint gcdr (const uint a, const uint b) { // greatest common divisor, Euclidean algorithm
	return b ? gcdr(b, a % b) : a; // recursion // Note: (uint) ? true{>0} : false{0};
}

uint gcd (uint a, uint b) {
	while (b) {
		uint t = a % b;
		a = b;
		b = t;
	}
	return a;
}

/* Notes on bitwise operations:
 * ((~i) & pow(2, k) - 1) === (!(i % pow(2, k)))
 * (i >> k) === (i / pow(2, k)) // right shift = division
 * (i << k) === (i * pow(2, k)) // left shift = multiplication
 * i >>= __builtin_ctz(i);  is equivalent to  while (!(i & 1)) i >>= 1;
 */

uint gcd2r (const uint a, const uint b) { // binary GCD algorithm or Stein's algorithm, the recursive version
	// simple cases, termination
	if (a == b || !a) // gcd(k, k) = k; gcd(0, b) = b
		return b;
	if (!b) // gcd(a, 0) = a
		return a;
	if (!(a & 1)) { // Case [1]: a is even
		if (!(b & 1)) // Case [1.1]: b is also even
			return gcd2r(a >> 1, b >> 1) << 1; // gcd(even, even) = 2 * gcd(even / 2, even / 2)
		else // Case [1.2]: b is odd
			return gcd2r(a >> 1, b); // [these are ...] gcd(even, odd) = gcd(even / 2, odd)
	} // Case [2]: a is odd
	if (!(b & 1)) // Case [2.1]: b is even
		return gcd2r(a, b >> 1); // [... the same] gcd(odd, even) = gcd(odd, even / 2)
	// Case [2.2]: both a and b are odd, gcd(odd, odd) = gcd((large - small) / 2, small)
	if (a > b)
		return gcd2r((a - b) >> 1, b);
	return gcd2r((b - a) >> 1, a); // reverse a and b, as a <= b
}

uint gcd2 (uint a, uint b) {
	if (!a)
		return b;
	if (!b)
		return a;
	int sh = __builtin_ctz(a | b); // not (a or b is odd) => (a and b are even)
	a >>= __builtin_ctz(a); // removed the greatest power of 2 which is a common factor of a and b
	do {
		b >>= __builtin_ctz(b); // b is now always odd
		if (a > b) { // make sure that b >= a, swap a and b if necessary
			uint t = a;
			a = b;
			b = t;
		}
		b -= a; // odd - odd = even; it is basically Euclidean algorithm but uses subtraction
	} while (b);
	return a << sh;
}

uint lcm (const uint a, const uint b) {
	return a / gcd(a, b) * b; // re-ordered to avoid overflowing even when the LCM does not
}

uint q005 (const uint n) {
	uint i, r = 1;
	for (i = 2; i <= n; i++)
		r = lcm(r, i);
	return r;
}

uint q005r (const uint n) {
	return n == 1 ? 1 : lcm(q005r(n - 1), n);
}

// question 5 by the sieve of Eratosthenes

uint* prime_esieve (const uint n) {
	uint* sieve = uintArrayBitCleared(n); // use uint (8-bit) as boolean (1-bit)
	uint i, j;
	for (i = 4; i <= n; i += 2)
		setBit(sieve, i); // set flag of each (even number > 2) to true
	for (i = 3; i <= sqrt(n); i += 2)
		if (!(getBit(sieve, i))) // i is prime
			for (j = i*i; j <= n; j += 2*i) // Skip i^2+i(k) = i(i+k) where k is odd since (i+k) is even
				setBit(sieve, j);
	return sieve;
}

uint q005_esieve (const uint n) { // n >= 2
	uint* sieve = prime_esieve(n);
	uint i, j = 2, p = 1; // j must be > 1
	for (i = 2; i <= n; i++) { // 1 is not sieved; find the nth primorial (product of the first n primes)
		if (getBit(sieve, i))
			continue;
		if (j == 1) {
			p *= i;
			continue;
		}
		j = log(n) / log(i);
		p *= pow(i, j);
	}
	free(sieve);
	return p;
}

// 6. square of sum minus minus sum of squares

uint sum_of_squares (const uint n) {
	return n * (n + 1) * (2*n + 1) / 6;
}

uint q006 (const uint n) {
	const uint sum_n = sum(n);
	return sum_n * sum_n - sum_of_squares(n);
}

// 7. the nth prime: given number of primes, find the value of prime
// 10. summation of primes: given (max value of primes + 1), find the sum of values

bool is_prime (const uint p) { // brute force trial division with some optimization tricks
	uint i;
	if (p <= 3) // 2, 3 are primes
		return (p > 1); // 1 is not prime; 0 is not a natural number
	if (!(p % 2) || !(p % 3)) // no even numbers nor multiples of 3
		return 0;
	for (i = 5; i <= sqrt(p); i += 6) // 6k +/- 1 from 5 to sqrt(p)
		if (!(p % i) || !(p % (i+2)))
			return 0;
	return 1;
}

uint q007 (const uint n) {
	if (n == 0)
		return 0; // invalid
	if (n <= 2)
		return n + 1; // primes: 2, 3, ...
	uint i = 3; // the 3rd prime is ...
	uint p = 5; // 5 (p is probable prime)
	bool q = 1; // 6k +/- 1 flag, set to true/+1 first [alternative: int q = 1;]
	while (i != n) {
		p += (q ? 2 : 4); // alternative: p += 3 - q;
		q = !q; // alternative: q = -q;
		if (is_prime(p))
			i++;
	}
	return p;
}

bool is_prime_cache (const uint p, const uint* primes) { // trial division with prime numbers stored
	uint i;
	if (p <= 1) // 1 is not prime; 0 is not a natural number
		return 0;
	// Removed {p = 2, returns true; even number, returns false} logic
	for (i = 1; primes[i] <= sqrt(p); i++)
		if (!(p % primes[i]))
			return 0;
	return 1;
}

uint q007_cache (const uint n) { // the first n primes
	const uint primesInit[] = { 3, 2, 3, 5 }; // element [0] is the counter, should also be the number of elements - 1

	if (n == 0)
		return 0; // invalid
	if (n <= 3)
		return primesInit[n]; // primes: 2, 3, 5 ...
	
	uint* primes = uintArray(n);
	memcpy(primes, primesInit, sizeof(primesInit));
	
	uint p = primes[primes[0]]; // 5 (p is probable prime)
	bool q = 1; // 6k +/- 1 flag, set to true/+1 first [alternative: int q = 1;]
	while (primes[0] != n) {
		p += (q ? 2 : 4); // alternative: p += 3 - q;
		q = !q; // alternative: q = -q;
		if (is_prime_cache(p, primes))
			primes[++primes[0]] = p; // also add p to cache
	}
	//return primes[primes[0]];
	free(primes);
	return p;
}

ulong q010_cache (const uint n) {
	const uint primesInit[] = { 3, 2, 3, 5 }; // element [0] is the counter, should also be the number of elements - 1

	uint i;
	ulong sum = 0;
	for (i = 1; i <= primesInit[0]; i++) {
		if (primesInit[i] >= n)
			return sum;
		sum += primesInit[i];
	}
	
	uint* primes = uintArray(n / log10(n)); // > (number of primes under n) > (n / log(n))
	memcpy(primes, primesInit, sizeof(primesInit));
	
	uint p = primes[primes[0]]; // 5 (p is probable prime)
	bool q = 1; // 6k +/- 1 flag, set to true/+1 first [alternative: int q = 1;]
	do {
		p += (q ? 2 : 4); // alternative: p += 3 - q;
		q = !q; // alternative: q = -q;
		if (p >= n)
			break;
		if (is_prime_cache(p, primes)) {
			primes[++primes[0]] = p; // also add p to cache
			sum += p;
		}
	} while (1);
	//return primes[0];
	//return primes[primes[0]];
	free(primes);
	return sum;
}

// question 7 and 10 by the sieve of Eratosthenes

uint* prime_esieve_odd (const uint n) { // 2x 8x memory usage improvement
	const uint size = (n - 1) / 2;
	uint* sieve = uintArrayBitCleared(size);
	uint i, j;
	for (i = 1; i <= (sqrt(n) - 1) / 2; i++) // index of sqrt(n)
		if (!getBit(sieve, i)) // (2i + 1) is prime
			for (j = 2*i*(i+1); j <= size; j += 2*i+1) // index of odd numbers squared = (2k+1)^2 = 2 [2i(i+1)] + 1
				setBit(sieve, j);
	return sieve;
}

ulong q010_esieve_odd (const uint n) { // n >= 3
	uint* sieve = prime_esieve_odd(n);
	ulong i = 0, p = 1, sum = 2; // 2 is the only even prime
	do {
		i++;
		p += 2;
		if (p >= n)
			break;
		if (!getBit(sieve, i))
			sum += p;
	} while (1);
	free(sieve);
	return sum;
}

uint q007_esieve_odd (const int n) { // n >= 2
	uint size;
	if (n <= 8601) // loser upper bound, error < +20%
		size = n * (log(n) + log(log(n)));
	else // tighter upper bound, error < +0.2%
		size = n * (log(n) + log(log(n)) - 0.9385);
	uint* sieve = prime_esieve_odd(size);
	uint i = 0, p = 1, c = 1; // 2 is the first, the only even prime
	do {
		i++;
		p += 2;
		if (!getBit(sieve, i))
			c++;
	} while (c != n);
	free(sieve);
	return p;
}

// 8. n adjacent digits in series having the greatest product

char* file_to_string (const char* path, const uint len) {
	FILE* fp;
	const uint lnlen = 80;
	char ln[lnlen], lns[lnlen];
	char* str = string(len);
	fp = fopen(path, "r");
	if (!fp)
		return str;
	while(fgets(ln, lnlen, fp) != NULL) {
		sscanf(ln, "%s\n", lns);
		strcat(str, lns);
	}
	fclose(fp);
	return str;
}

uint q008 (const char* path, const uint len, const uint d) {
	uint i, j;
	uint max = 0;
	char* str = file_to_string(path, len);
	for (i = 0; i <= len - d; i++) {
		uint p;
		p = 1;
		for (j = 0; j < d; j++)
			p *= charToInt(str[i + j]);
		if (p > max)
			max = p;
	}
	return max;
}

ulong q008_fast (const char* path, const uint len, const uint d) { // dynamic programming
	uint i, z = 0; // number of zeros in product
	ulong p = 1, max = 0;
	char* str = file_to_string(path, len);
	for (i = 0; i < len; i++) {
		// printf("p=%ld, max=%ld, z=%u\n", p, max, z);
		if (i+1 > d) { // i [0 .. d - 1]
			if (str[i - d] == '0')
				z--;
			else
				p /= charToInt(str[i - d]);
		}
		if (str[i] == '0')
			z++;
		else
			p *= charToInt(str[i]);
		if (i+1 >= d && !z && p > max)
			max = p;
	}
	return max;
}

// 9. Pythagorean triplet given perimeter

uint q009 (const uint s) {
	uint a, b, c;
	for (a = 3; a <= (s - 3) / 3; a++) { // a + (a+1) + [(a+1) + 1] <= s
		for (b = a + 1; b <= (s - a - 1) / 2; b++) { // a + b + (b+1) <= s
			c = s - a - b;
			if (a*a + b*b == c*c)
				return a * b * c;
		}
	}
	return 0;
}

uint q009_fast (const uint s) { // s is perimeter which must be even
	const uint s2 = s / 2; // s = a+b+c = d(2m^2 + 2mn) = 2(d)(m)(m+n) where (n > 0), (m + n > m)
	uint m, s2m, k, d, n, a, b, c;
	for (m = 2; m <= sqrt(s); m++) {
		if (!(s2 % m)) { // m is divisible by s/2
			s2m = s2 / m;
			s2m >>= __builtin_ctz(s2m); // remove all prime factors 2 (binary: trailing zeros)
			if (m & 1) // m is odd
				k = m + 2; // n is even, k = odd + even = odd
			else // m is even
				k = m + 1; // n is odd, k = even + odd = odd
			while (k < 2*m && k <= s2m) { // (n < m) and (1 <= d)
				if (!(s2m % k) && gcd(k, m) == 1) { // k is divisible by s/2/m; (m + n) and m coprime
					d = s2m / k;
					n = k - m;
					a = d * (m*m - n*n);
					b = d * 2 * m * n;
					c = d * (m*m + n*n);
					return a * b * c;
				}
				k += 2;
			}
		}
	}
	return 0;
}

// for question 11 later
uchar* file_to_uchar_2d (const char* path, const uint side) { // 2-digit uint (uchar) separated by 1 character
	FILE* fp;
	const uint lnlen = 3*side + 1;
	char ln[lnlen];
	uchar* a = ucharArray(side*side);
	fp = fopen(path, "r");
	if (!fp)
		return a;
	uint i, ios = 0;
	while(fgets(ln, lnlen, fp) != NULL) {
		for (i = 0; i < side; i++) {
			a[ios + i] = 10 * charToInt(ln[3*i]) + charToInt(ln[3*i + 1]);
		}
		ios += side;
	}
	fclose(fp);
	return a;
}

// for question 16 later
ulong sum_of_digits (ulong n) {
	ulong sum = 0;
	while (n) {
		sum += n % 10;
		n /= 10;
	}
	return sum;
}

ulong factorial (const uint n) { // unsigned long int overflows when n > 25
	return n == 1 ? 1 : n * factorial(n - 1);
}

int main (int args, char* argv[]) {
	const long result[] = {
		// 1. sum of multiples of a or b inclusive
		q001_fast(3, 5, 1000),
		// 2. sum of even Fibonacci numbers
		q002_fast(4 * pow(10,6)),
		// 3. largest prime factor (prime sieve)
		q003(600851475143),
		// 4. largest (even-digit) palindrome product of two n-digit numbers
		q004_fast(3),
		// 5. LCM of 1 to n
		q005_esieve(20),
		// 6. square of sum minus sum of squares
		q006(100),
		// 7. the nth prime: given number of primes, find the value of prime
		q007_esieve_odd(10001),
		// 8. n adjacent digits in series having the greatest product
		q008_fast("input/q008.txt", 1000, 13),
		// 9. Pythagorean triplet given perimeter
		q009_fast(1000),
		// 10. summation of primes: given (max value of primes + 1), find the sum of values
		q010_esieve_odd(2 * pow(10,6)),
	};
	uint i, len = sizeof(result) / sizeof(long);
	for (i = 0; i < len; i++)
		printf("%u. %ld \n", i+1, result[i]);
	
	for (i = 0; i < args; i++) {
		printf("%s ", argv[i]);
	}
	printf("\n");
	return 0;
}
