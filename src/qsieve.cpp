/*
Copyright 2022, Yves Gallot

qsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <queue>
#include <thread>
#include <mutex>

static uint64_t neg_mod(const uint64_t x, const uint64_t p) { return (x == 0) ? 0 : p - x; }

static uint64_t add_mod(const uint64_t x, const uint64_t y, const uint64_t p)
{
	const uint64_t c = (x >= p - y) ? p : 0;
	return x + y - c;
}

static uint64_t sub_mod(const uint64_t x, const uint64_t y, const uint64_t p)
{
	const uint64_t c = (x < y) ? p : 0;
	return x - y + c;
}

static uint64_t half_mod(const uint64_t x, const uint64_t p)
{
	const uint64_t c = (x % 2 != 0) ? p : 0;
	return (x + c) / 2;
}

class mod
{
private:
	const uint64_t _p;
	const int _shift;
	const uint64_t _q;

private:
	constexpr static int log2_64(const uint64_t x) { return 63 - __builtin_clzll(x); }
	constexpr static uint64_t mulhi_64(const uint64_t x, const uint64_t y) { return uint64_t((x * __uint128_t(y)) >> 64); }

public:
	mod(const uint64_t p) : _p(p), _shift(log2_64(p) - 1), _q(uint64_t((__uint128_t(1) << (64 + _shift)) / p)) {}

	// 2 < p
	constexpr uint64_t oneHalf() const { return (_p + 1) / 2; }

	// 3 < p < 2^63
	constexpr uint64_t oneThird() const
	{
		return (_p % 3 == 2) ? (_p + 1) / 3 : (2 * _p + 1) / 3;
	}

	// 5 < p < 2^64 / 3
	constexpr uint64_t oneFifth() const
	{
		const uint64_t r = _p % 5;
		if (r == 1) return (4 * _p + 1) / 5;
		if (r == 2) return (2 * _p + 1) / 5;
		if (r == 3) return (3 * _p + 1) / 5;
		return (_p + 1) / 5;
	}

	uint64_t neg(const uint64_t x) const { return neg_mod(x, _p); }
	uint64_t add(const uint64_t x, const uint64_t y) const { return add_mod(x, y, _p); }
	uint64_t sub(const uint64_t x, const uint64_t y) const { return sub_mod(x, y, _p); }
	uint64_t half(const uint64_t x) const { return half_mod(x, _p); }

	// Barrett's product: let n = 63, r = ceil(log2(p)) [2^{r - 1} < p < 2^r], s = r - 2 = floor(log2(p)) - 1,
	// t = n + 1 = 64, q = floor(2^{s + t} / p). Then the number of iterations h = 1.
	// We must have x.y < alpha.p with alpha = 2^{n-2}. If p <= 2^{n-2} = 2^61 then x^2 < p^2 <= alpha.p.
	uint64_t mul(const uint64_t x, const uint64_t y) const
	{
		const __uint128_t xy = x * __uint128_t(y);		// 0 <= ab < 2^{2r}
		const uint64_t q_p = uint64_t(xy >> _shift);	// 0 <= q_p < 2^{r + 2} <= 2^63
		const uint64_t r = uint64_t(xy - mulhi_64(q_p, _q) * __uint128_t(_p));
		return (r >= _p) ? r - _p : r;
	}

	uint64_t oneHalfPow(const uint64_t e) const
	{
		uint64_t r = oneHalf();
		for (int b = log2_64(e) - 1; b >= 0; --b)
		{
			r = mul(r, r);
			if ((e & (uint64_t(1) << b)) != 0) r = half(r);
		}
		return r;
	}

	uint64_t twoPow(const uint64_t e) const
	{
		uint64_t r = 2;
		for (int b = log2_64(e) - 1; b >= 0; --b)
		{
			r = mul(r, r);
			if ((e & (uint64_t(1) << b)) != 0) r = add(r, r);
		}
		return r;
	}

	bool spsp2() const
	{
		// n - 1 = 2^k * r
		uint64_t r = _p - 1;
		int k = 0; for (; r % 2 == 0; r /= 2) ++k;

		uint64_t x = twoPow(r);
		if (x == 1) return true;

		// Compute x^(2^i) for 0 <= i < n.  If any are -1, n is a a-spsp.
		for (; k > 0; --k)
		{
			if (x == _p - 1) return true;
			x = mul(x, x);
		}

		return false;
	}
};

class bitmap
{
private:
	const size_t _k_size, _n_size;
	std::vector<bool> _bmp;

public:
	bitmap(const size_t k_size, const size_t n_size) : _k_size(k_size), _n_size(n_size)
	{
		_bmp.resize(k_size * n_size, false);
	}

	bool get(const size_t k, const size_t n) const { return _bmp[n * _k_size + k]; }
	void set(const size_t k, const size_t n) { if (!_bmp[n * _k_size + k]) _bmp[n * _k_size + k] = true; }

	size_t size() const { return _bmp.size(); }

	size_t count() const
	{
		size_t cnt = 0;
		for (bool b : _bmp) cnt += b ? 0 : 1;
		return cnt;
	}
};

class qsieve
{
private:
	// From 3321925 to 3371925: zero-padded FMA3 FFT length 336K
	const uint64_t _n_min = 3321925;
	const size_t _n_range = 64;
	const uint64_t _p_min, _p_max, _k_min, _k_max;
	const size_t _k_range;

	static const size_t p_size = 1024;

	struct PArray
	{
		uint64_t p[p_size];	// 8 KB
	};

	struct PR { uint64_t p, r; };

	struct PRArray
	{
		PR pr[p_size];
	};

	static const size_t max_queue_size = 1024;

	std::mutex _p_queue_mutex;
	std::queue<PArray> _p_queue;

	std::mutex _pr_queue_mutex;
	std::queue<PRArray> _pr_queue;

	bool _end_p = false, _end_r = false;

private:
	void gen_p()
	{
		// Segmented sieve of Eratosthenes: outputs have no factor < 65537.
		static const uint32_t sp_max = 1 << 16;
		static const size_t sieve_size = sp_max / 2;	// sieve with an odd prime table.
		static const size_t odd_prime_count = 6541;		// # odd primes with p < 2^16.

		bool sieve[sieve_size];
		uint32_t prm[odd_prime_count];
		uint32_t prm_ptr[odd_prime_count];

		prm[0] = 3; prm[1] = 5; prm[2] = 7;
		uint32_t i = 3;
		for (uint32_t k = 11; k < sp_max; k += 2)
		{
			const uint32_t s = uint32_t(std::sqrt(double(k))) + 1;
			uint32_t d; for (d = 3; d <= s; d += 2) if (k % d == 0) break;
			if (d > s) prm[i++] = k;
		}

		// if (i != odd_prime_count) throw;

		for (size_t k = 0; k < sieve_size; ++k) sieve[k] = false;

		const uint64_t p0 = (_p_min / sp_max) * sp_max, p1 = (_p_max / sp_max + 1) * sp_max;
		std::cout << "p in [" << p0 << "; " << p1 << "] " << std::endl;

		if (p0 == 0)
		{
			sieve[0] = true;	// p = 1
			for (size_t i = 0; i < odd_prime_count; ++i)
			{
				const size_t p = prm[i];
				bool first = true;
				size_t k = (p - 1) / 2;
				for (; k < sieve_size; k += p) if (first) first = false; else sieve[k] = true;
				prm_ptr[i] = uint32_t(k - sieve_size);
			}
		}
		else
		{
			for (size_t i = 0; i < odd_prime_count; ++i)
			{
				const size_t p = prm[i];
				size_t o = p - size_t(p0 % p); if (o % 2 == 0) o += p;
				size_t k = (o - 1) / 2;
				for (; k < sieve_size; k += p) sieve[k] = true;
				prm_ptr[i] = uint32_t(k - sieve_size);
			}
		}

		PArray p_array;
		size_t p_array_i = 0;
		size_t queue_size = 0;

		for (uint64_t jp = p0; true; jp += sp_max)
		{
			for (size_t kp = 0; kp < sieve_size; ++kp)
			{
				if (!sieve[kp])
				{
					const uint64_t p = jp + 2 * kp + 1;
					if ((p > 5) && ((p < 65537 * uint64_t(65537)) || mod(p).spsp2()))
					{
						p_array.p[p_array_i] = p;
						p_array_i = (p_array_i + 1) % p_size;
						if (p_array_i == 0)
						{
							std::lock_guard<std::mutex> guard(_p_queue_mutex);
							_p_queue.push(p_array);
							queue_size = _p_queue.size();
							if (p >= p1)
							{
								_end_p = true;
								return;
							}
						}
					}
				}
			}

			for (size_t k = 0; k < sieve_size; ++k) sieve[k] = false;

			for (size_t i = 0; i < odd_prime_count; ++i)
			{
				size_t k = prm_ptr[i], p = prm[i];
				for (; k < sieve_size; k += p) sieve[k] = true;
				prm_ptr[i] = uint32_t(k - sieve_size);
			}

			while (queue_size > max_queue_size)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(100));
				std::lock_guard<std::mutex> guard(_p_queue_mutex);
				queue_size = _p_queue.size();
			}
		}
	}

	void gen_r()
	{
		PArray p_array;
		PRArray pr_array;
		size_t queue_size = 0;

		while (true)
		{
			bool found = false;
			{
				std::lock_guard<std::mutex> guard(_p_queue_mutex);
				if (!_p_queue.empty())
				{
					found = true;
					p_array = _p_queue.front();
					_p_queue.pop();
				}
			}

			if (!found)
			{
				if (_end_p) { _end_r = true; return; }
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
			}
			else
			{
				for (size_t i = 0; i < p_size; ++i)
				{
					const uint64_t p = p_array.p[i];
					const mod m = mod(p);
					PR & pr = pr_array.pr[i];
					pr.p = p;
					pr.r = m.mul(m.oneHalfPow(_n_min - 1), m.mul(m.oneThird(), m.oneFifth()));	// r = 1/2^(n_min - 1) / 15 (mod p)
				}

				std::lock_guard<std::mutex> guard(_pr_queue_mutex);
				_pr_queue.push(pr_array);
				queue_size = _pr_queue.size();
			}

			while (queue_size > max_queue_size)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(100));
				std::lock_guard<std::mutex> guard(_pr_queue_mutex);
				queue_size = _pr_queue.size();
			}
		}
	}

	static uint64_t kpos(const uint64_t p, const uint64_t k_0, const uint64_t r) { return half_mod(sub_mod(r, k_0, p), p); }
	static uint64_t kneg(const uint64_t p, const uint64_t k_0, const uint64_t r) { return half_mod(neg_mod(add_mod(r, k_0, p), p), p); }

	void sieve(bitmap & bmap, const uint64_t p, const uint64_t r_0)
	{
		const size_t k_range = _k_range;
		const uint64_t k_0 = (_k_min / 15) % p;

		uint64_t r = r_0;	// r = 1/2^(n_min - 1) / 15 (mod p)

		for (uint64_t k = kpos(p, k_0, r); k < k_range; k += p)
		{
			bmap.set(k, -1 + 1);
		}
		r = half_mod(r, p);	// r = 1/2^n_min / 15 (mod p)

		for (uint64_t k = kpos(p, k_0, r); k < k_range; k += p)
		{
			bmap.set(k, 0);
			bmap.set(k, 0 + 1);
		}
		for (uint64_t k = kneg(p, k_0, r); k < k_range; k += p)
		{
			bmap.set(k, 0);
		}
		r = half_mod(r, p);	// 1/2^(n_min + 1) / 15 (mod p)

		for (size_t n = 1; n < _n_range - 1; ++n)
		{
			for (uint64_t k = kpos(p, k_0, r); k < k_range; k += p)
			{
				bmap.set(k, n - 1);
				bmap.set(k, n);
				bmap.set(k, n + 1);
			}
			for (uint64_t k = kneg(p, k_0, r); k < k_range; k += p)
			{
				bmap.set(k, n);
			}
			r = half_mod(r, p);	// 1/2^(n_min + n) / 15 (mod p)
		}

		for (uint64_t k = kpos(p, k_0, r); k < k_range; k += p)
		{
			bmap.set(k, _n_range - 1 - 1);
			bmap.set(k, _n_range - 1);
		}
		for (uint64_t k = kneg(p, k_0, r); k < k_range; k += p)
		{
			bmap.set(k, _n_range - 1);
		}
		r = half_mod(r, p);	// 1/2^(n_min + _n_range) / 15 (mod p)

		for (uint64_t k = kpos(p, k_0, r); k < k_range; k += p)
		{
			bmap.set(k, _n_range - 1);
		}
	}

public:
	qsieve(const uint64_t p_max, const uint64_t k_min, const uint64_t k_max) : _p_min(2), _p_max(p_max), _k_min(k_min), _k_max(k_max),
		_k_range(size_t(k_max - k_min) / 30 + 1)
	{
		std::cout << "k_min = " << k_min << ", k_max = " << k_max << ", n_min = " << _n_min << ", n_max = " << _n_min + _n_range - 1 << std::endl;

		std::thread t_gen_p([=] { gen_p(); }); t_gen_p.detach();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		std::thread t_gen_r([=] { gen_r(); }); t_gen_r.detach();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));

		bitmap bmap(_k_range, _n_range);

		std::cout << bmap.size() / (8 << 20) << " MB" << std::endl;

		PRArray pr_array;
		uint64_t last_p = 0;
		double duration = 0;
		auto t0 = std::chrono::steady_clock::now();

		while (true)
		{
			bool found = false;
			{
				std::lock_guard<std::mutex> guard(_pr_queue_mutex);
				if (!_pr_queue.empty())
				{
					found = true;
					pr_array = _pr_queue.front();
					_pr_queue.pop();
				}
			}

			if (!found)
			{
				if (_end_r) break;
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
			}
			else
			{
				for (size_t i = 0; i < p_size; ++i)
				{
					const PR & pr = pr_array.pr[i];
					const uint64_t p = pr.p, r = pr.r;
					if (duration == 0) { std::cout << p << "\r"; std::cout.flush(); }
					if (p < p_max) { last_p = p; sieve(bmap, p, r); }
				}

				const auto t1 = std::chrono::steady_clock::now();
				const std::chrono::duration<double> dt = t1 - t0;
				if ((dt.count() > 5 * 60) || (duration == 0))
				{
					t0 = t1; duration += dt.count();
					const size_t expected = size_t(0.41252 * (k_max - k_min + 1) * _n_range / std::pow(std::log(double(last_p)), 4));
					std::cout << std::scientific << std::setprecision(2) << double(last_p) << ": ";
					std::cout << bmap.count() << " candidates, " << expected << " expected, " << std::lrint(duration) << " sec." << std::endl;
				}
			}
		}

		const std::chrono::duration<double> dt = std::chrono::steady_clock::now() - t0;
		duration += dt.count();
		const size_t expected = size_t(0.41252 * (k_max - k_min + 1) * _n_range / std::pow(std::log(double(p_max)), 4));
		std::cout << std::scientific << std::setprecision(2) << double(p_max) << ": ";
		std::cout << bmap.count() << " candidates, " << expected << " expected, " << std::lrint(duration) << " sec." << std::endl;

		for (size_t k = 0; k < _k_range; ++k)
		{
			for (size_t n = 0; n < _n_range; ++n)
			{
				if (!bmap.get(k, n))
				{
					std::cout << k_min + 30 * k << ", " << _n_min + n << std::endl;
				}
			}
		}
	}

	virtual ~qsieve() {}
};

int main(int argc, char * argv[])
{
	std::cerr << "qsieve: quad sieve for twin and Sophie Germain primes" << std::endl;
	std::cerr << " Copyright (c) 2022, Yves Gallot" << std::endl;
	std::cerr << " qsieve is free source code, under the MIT license." << std::endl << std::endl;
	std::cerr << " Usage: qsieve <p_max> <k_min> <k_max>" << std::endl << std::endl;

	uint64_t p_max = (argc > 1) ? std::atoll(argv[1]) : uint64_t(-1) / 4;
	if (p_max < 7) p_max = 7;
	if (p_max > uint64_t(-1) / 4) p_max = uint64_t(-1) / 4;

	uint64_t k_min = (argc > 2) ? std::atoll(argv[2]) : 15;
	k_min /= 30; k_min *= 30; k_min += 15;

	uint64_t k_max = (argc > 3) ? std::atoll(argv[3]) : k_min + 30 * uint64_t(143165576);
	k_max /= 30; k_max *= 30; k_max += 15;
	if (k_max < k_min) k_max = k_min;

	qsieve(p_max, k_min, k_max);

	return EXIT_SUCCESS;
}
