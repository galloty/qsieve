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

constexpr static uint64_t neg_mod(const uint64_t x, const uint64_t p) { return (x == 0) ? 0 : p - x; }

constexpr static uint64_t add_mod(const uint64_t x, const uint64_t y, const uint64_t p)
{
	const uint64_t c = (x >= p - y) ? p : 0;
	return x + y - c;
}

constexpr static uint64_t sub_mod(const uint64_t x, const uint64_t y, const uint64_t p)
{
	const uint64_t c = (x < y) ? p : 0;
	return x - y + c;
}

constexpr static uint64_t half_mod(const uint64_t x, const uint64_t p)
{
	const uint64_t c = (x % 2 != 0) ? p : 0;
	return (x + c) / 2;
}

// 3 < p < 2^63
constexpr static uint64_t oneThird_mod(const uint64_t p)
{
	return (p % 3 == 2) ? (p + 1) / 3 : (2 * p + 1) / 3;
}

// 5 < p < 2^64 / 3
constexpr static uint64_t oneFifth_mod(const uint64_t p)
{
	const uint64_t r = p % 5;
	if (r == 1) return (4 * p + 1) / 5;
	if (r == 2) return (2 * p + 1) / 5;
	if (r == 3) return (3 * p + 1) / 5;
	return (p + 1) / 5;
}

class res4
{
private:
	uint64_t _r[4];

public:
	uint64_t operator [](const size_t i) const { return _r[i]; }
	uint64_t & operator [](const size_t i) { return _r[i]; }
};

class mod4
{
private:
	res4 _p;
	int _shift[4];
	uint64_t _q[4];

private:
	constexpr static int log2_64(const uint64_t x) { return 63 - __builtin_clzll(x); }
	constexpr static uint64_t mulhi_64(const uint64_t x, const uint64_t y) { return uint64_t((x * __uint128_t(y)) >> 64); }

	res4 half(const res4 & x) const
	{
		res4 r;
		for (size_t i = 0; i < 4; ++i) r[i] = half_mod(x[i], _p[i]);
		return r;
	}

	// Barrett's product: let n = 63, r = ceil(log2(p)) [2^{r - 1} < p < 2^r], s = r - 2 = floor(log2(p)) - 1,
	// t = n + 1 = 64, q = floor(2^{s + t} / p). Then the number of iterations h = 1.
	// We must have x.y < alpha.p with alpha = 2^{n-2}. If p <= 2^{n-2} = 2^61 then x^2 < p^2 <= alpha.p.
	uint64_t mul(const uint64_t x, const uint64_t y, const size_t i) const
	{
		const __uint128_t xy = x * __uint128_t(y);		// 0 <= ab < 2^{2r}
		const uint64_t q_p = uint64_t(xy >> _shift[i]);	// 0 <= q_p < 2^{r + 2} <= 2^63
		const uint64_t r = uint64_t(xy - mulhi_64(q_p, _q[i]) * __uint128_t(_p[i]));
		return (r >= _p[i]) ? r - _p[i] : r;
	}

public:
	mod4(const res4 & p)
	{
		for (size_t i = 0; i < 4; ++i)
		{
			_p[i] = p[i];
			_shift[i] = log2_64(p[i]) - 1;
			_q[i] = uint64_t((__uint128_t(1) << (64 + _shift[i])) / p[i]);
		}
	}

	res4 mul(const res4 & x, const res4 & y) const
	{
		res4 r;
		for (size_t i = 0; i < 4; ++i) r[i] = mul(x[i], y[i], i);
		return r;
	}

	res4 oneFifteenth() const
	{
		res4 r;
		for (size_t i = 0; i < 4; ++i) r[i] = mul(oneThird_mod(_p[i]), oneFifth_mod(_p[i]), i);
		return r;
	}

	res4 oneHalf_pow(const uint64_t e) const
	{
		res4 r;
		for (size_t i = 0; i < 4; ++i) r[i] = (_p[i] + 1) / 2;
		int b = log2_64(e) - 1;
		do
		{
			r = mul(r, r);
			if ((e & (uint64_t(1) << b)) != 0) r = half(r);
		} while (--b >= 0);
		return r;
	}
};

class bitmap64
{
private:
	const size_t _k_size;
	uint64_t * const _bmp;

public:
	bitmap64(const size_t k_size) : _k_size(k_size), _bmp(new uint64_t[k_size])
	{
		for (size_t k = 0; k < k_size; ++k) _bmp[k] = 0;
	}
	virtual ~bitmap64() { delete[] _bmp; }

	bool get(const size_t k, const size_t n) const { const uint64_t m = uint64_t(1) << n; return (_bmp[k] & m) != 0; }

	void set(const size_t k, const size_t n) { const uint64_t m = uint64_t(1) << n; if ((_bmp[k] & m) != m) _bmp[k] |= m; }
	void set2(const size_t k, const size_t n) { const uint64_t m = uint64_t(3) << n; if ((_bmp[k] & m) != m) _bmp[k] |= m; }
	void set3(const size_t k, const size_t n) { const uint64_t m = uint64_t(7) << n; if ((_bmp[k] & m) != m) _bmp[k] |= m; }

	size_t size() const { return _k_size * 8 * sizeof(uint64_t); }

	size_t count() const
	{
		size_t cnt = size();
		for (size_t k = 0, k_size = _k_size; k < k_size; ++k) cnt -= __builtin_popcountll(_bmp[k]);
		return cnt;
	}
};

template<typename T>
class fifo
{
private:
	static const size_t max_queue_size = 1024;

	std::mutex _mutex;
	std::queue<T> _queue;
	bool _end = false;

public:
	void end() { std::lock_guard<std::mutex> guard(_mutex); _end = true; }

	void push(const T & val)
	{
		_mutex.lock();
		while (_queue.size() >= max_queue_size)
		{
			_mutex.unlock();
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
			_mutex.lock();
		}
		_queue.push(val);
		_mutex.unlock();
	}

	bool pop(T & val)
	{
		_mutex.lock();
		while (_queue.empty())
		{
			if (_end)
			{
				_mutex.unlock();
				return false;
			}
			_mutex.unlock();
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
			_mutex.lock();
		}
		val = _queue.front();
		_queue.pop();
		_mutex.unlock();
		return true;
	}
};

class qsieve
{
private:
	// From 3321925 to 3371925: zero-padded FMA3 FFT length 336K
	const uint64_t _n_min = 3321925;
	const size_t _n_range = 64;
	const uint64_t _p_min, _p_max, _k_min_15;
	const size_t _k_range;

	static const size_t p_size = 1024;

	struct PArray { uint64_t p[p_size]; };

	struct PR { uint64_t p, r; };
	struct PRArray { PR pr[p_size]; };

	fifo<PArray> _p_queue;
	fifo<PRArray> _pr_queue;

private:
	void gen_p()
	{
		// Segmented sieve of Eratosthenes: outputs have no factor < 2^20. 2^40 ~ 10^12
		static const uint32_t sp_max = 1 << 20;
		static const size_t sieve_size = sp_max / 2;	// sieve with an odd prime table.

		std::vector<bool> sieve(sieve_size, false);
		std::vector<uint32_t> v_prm, v_prm_ptr;

		v_prm.push_back(3); v_prm.push_back(5); v_prm.push_back(7);
		for (uint32_t k = 11; k < sp_max; k += 2)
		{
			const uint32_t s = uint32_t(std::sqrt(double(k))) + 1;
			uint32_t d; for (d = 3; d <= s; d += 2) if (k % d == 0) break;
			if (d > s) v_prm.push_back(k);
		}

		const uint64_t p0 = (_p_min / sp_max) * sp_max, p1 = (_p_max / sp_max + 1) * sp_max;
		std::cout << "Sieve of Eratosthenes: " << v_prm.size() << " primes (" << v_prm.size() * 8 / 1024 << " kB), bitmap size: " << sieve_size / (8 * 1024)
				  << " kB, p in [" << std::max(p0, uint64_t(7)) << "; " << p1 << "] " << std::endl;

		if (p0 == 0)
		{
			sieve[0] = true;	// p = 1
			for (const size_t p : v_prm)
			{
				bool first = true;
				size_t k = (p - 1) / 2;
				for (; k < sieve_size; k += p) if (first) first = false; else sieve[k] = true;
				v_prm_ptr.push_back(uint32_t(k - sieve_size));
			}
		}
		else
		{
			for (const size_t p : v_prm)
			{
				size_t o = p - size_t(p0 % p); if (o % 2 == 0) o += p;
				size_t k = (o - 1) / 2;
				for (; k < sieve_size; k += p) sieve[k] = true;
				v_prm_ptr.push_back(uint32_t(k - sieve_size));
			}
		}

		const size_t odd_prime_count = v_prm.size();
		const uint32_t * const prm = v_prm.data();
		uint32_t * const prm_ptr = v_prm_ptr.data();

		PArray p_array;
		size_t p_array_i = 0;

		for (uint64_t jp = p0; true; jp += sp_max)
		{
			for (size_t kp = 0; kp < sieve_size; ++kp)
			{
				if (!sieve[kp])
				{
					const uint64_t p = jp + 2 * kp + 1;
					if (p > 5)
					{
						p_array.p[p_array_i] = p;
						p_array_i = (p_array_i + 1) % p_size;
						if (p_array_i == 0)
						{
							_p_queue.push(p_array);
							if (p >= p1)
							{
								_p_queue.end();
								return;
							}
						}
					}
				}
			}

			sieve.assign(sieve_size, false);

			for (size_t i = 0; i < odd_prime_count; ++i)
			{
				size_t k = prm_ptr[i], p = prm[i];
				for (; k < sieve_size; k += p) sieve[k] = true;
				prm_ptr[i] = uint32_t(k - sieve_size);
			}
		}
	}

	void gen_r()
	{
		while (true)
		{
			PArray p_array;
			if (!_p_queue.pop(p_array)) { _pr_queue.end(); break; }

			PRArray pr_array;
			for (size_t i = 0; i < p_size; i += 4)
			{
				res4 p; for (size_t j = 0; j < 4; ++j) p[j] = p_array.p[i + j];
				const mod4 m = mod4(p);

				// r = 1/2^(n_min - 1) / 15 (mod p)
				const res4 r = m.mul(m.oneHalf_pow(_n_min - 1), m.oneFifteenth());
				for (size_t j = 0; j < 4; ++j)
				{
					PR & pr = pr_array.pr[i + j];
					pr.p = p[j]; pr.r = r[j];
				}
			}
			_pr_queue.push(pr_array);
		}
	}

	static size_t kpos(const uint64_t p, const uint64_t k_0, const uint64_t r) { return size_t(half_mod(sub_mod(r, k_0, p), p)); }
	static size_t kneg(const uint64_t p, const uint64_t k_0, const uint64_t r) { return size_t(half_mod(neg_mod(add_mod(r, k_0, p), p), p)); }

	void fill_sieve(bitmap64 & bmap, const uint64_t p, const uint64_t r_0)
	{
		const uint64_t k_0 = _k_min_15 % p;
		const size_t k_range = _k_range;

		uint64_t r = r_0;	// r = 1/2^(n_min - 1) / 15 (mod p)

		for (size_t k = kpos(p, k_0, r); k < k_range; k += p) bmap.set(k, 0); // -1 - 1, -1, -1 + 1
		// for (size_t k = kneg(p, k_0, r); k < k_range; k += p); // -1
		r = half_mod(r, p);	// r = 1/2^n_min / 15 (mod p)

		for (size_t k = kpos(p, k_0, r); k < k_range; k += p) bmap.set2(k, 0); // 0 - 1, 0, 0 + 1
		for (size_t k = kneg(p, k_0, r); k < k_range; k += p) bmap.set(k, 0); // 0
		r = half_mod(r, p);	// 1/2^(n_min + 1) / 15 (mod p)

		for (size_t n = 1; n < _n_range - 1; ++n)
		{
			for (size_t k = kpos(p, k_0, r); k < k_range; k += p) bmap.set3(k, n - 1); // n - 1, n, n + 1
			for (size_t k = kneg(p, k_0, r); k < k_range; k += p) bmap.set(k, n); // n
			r = half_mod(r, p);	// 1/2^(n_min + n) / 15 (mod p)
		}

		for (size_t k = kpos(p, k_0, r); k < k_range; k += p) bmap.set2(k, _n_range - 2); // (_n_range - 1) - 1, (_n_range - 1), (_n_range - 1) + 1
		for (size_t k = kneg(p, k_0, r); k < k_range; k += p) bmap.set(k, _n_range - 1); // (_n_range - 1)
		r = half_mod(r, p);	// 1/2^(n_min + n_range) / 15 (mod p)

		for (size_t k = kpos(p, k_0, r); k < k_range; k += p) bmap.set(k, _n_range - 1); // _n_range - 1, _n_range, _n_range + 1
		// for (size_t k = kneg(p, k_0, r); k < k_range; k += p); // _n_range
	}

public:
	qsieve(const uint64_t p_max, const uint64_t k_min, const uint64_t k_max) : _p_min(2), _p_max(p_max), _k_min_15(k_min / 15), _k_range(size_t(k_max - k_min) / 30 + 1)
	{
		std::cout << "k_min = " << k_min << ", k_max = " << k_max << ", n_min = " << _n_min << ", n_max = " << _n_min + _n_range - 1 << std::endl;

		std::thread t_gen_p([=] { gen_p(); }); t_gen_p.detach();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		std::thread t_gen_r([=] { gen_r(); }); t_gen_r.detach();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));

		bitmap64 bmap(_k_range);	// _n_range must be 64

		std::cout << "Bitmap size: " << bmap.size() / (8 << 20) << " MB" << std::endl;

		uint64_t last_p = 0, last_p0 = 0;
		double duration = 0;
		auto t0 = std::chrono::steady_clock::now();

		while (true)
		{
			PRArray pr_array;
			if (!_pr_queue.pop(pr_array)) break;

			for (size_t i = 0; i < p_size; ++i)
			{
				const PR & pr = pr_array.pr[i];
				const uint64_t p = pr.p, r = pr.r;
				if (duration == 0) { std::cout << p << "\r"; std::cout.flush(); }
				if (p < p_max) { last_p = p; fill_sieve(bmap, p, r); }
			}

			const auto t1 = std::chrono::steady_clock::now();
			const double dt = std::chrono::duration<double>(t1 - t0).count();
			if ((dt > 15 * 60) || (duration == 0))
			{
				std::cout << std::scientific << std::setprecision(2) << double(last_p) << " (+" << double(last_p - last_p0) << "): ";
				duration += dt; last_p0 = last_p; t0 = t1;
				const size_t expected = size_t(0.41252 * (k_max - k_min + 1) * _n_range / std::pow(std::log(double(last_p)), 4));
				std::cout << bmap.count() << " candidates, " << expected << " expected, " << std::lrint(duration) << " sec." << std::endl;
			}
		}

		const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
		std::cout << std::scientific << std::setprecision(2) << double(p_max) << " (+" << double(p_max - last_p0) << "): ";
		duration += dt;
		const size_t expected = size_t(0.41252 * (k_max - k_min + 1) * _n_range / std::pow(std::log(double(p_max)), 4));
		std::cout << bmap.count() << " candidates, " << expected << " expected, " << std::lrint(duration) << " sec." << std::endl;

		std::ofstream outFile("qsieve.log");
		for (size_t k = 0, k_range = _k_range; k < k_range; ++k)
		{
			for (size_t n = 0; n < _n_range; ++n)
			{
				if (!bmap.get(k, n))
				{
					outFile << k_min + 30 * k << "\t" << _n_min + n << std::endl;
				}
			}
		}
		outFile.close();
	}

	virtual ~qsieve() {}
};

int main(int argc, char * argv[])
{
	std::cerr << "qsieve: quad sieve for twin and Sophie Germain primes" << std::endl;
	std::cerr << " Copyright (c) 2022, Yves Gallot" << std::endl;
	std::cerr << " qsieve is free source code, under the MIT license." << std::endl << std::endl;
	std::cerr << " Usage: qsieve <p_max> <k_min> <k_max>" << std::endl << std::endl;

	// uint64_t p_max = (argc > 1) ? std::atoll(argv[1]) : uint64_t(-1) / 4;
	// uint64_t T = 1; for (size_t i = 0; i < 8; ++i) T *= 10;	// 6705 candidates for 1000000
	uint64_t T = 1; for (size_t i = 0; i < 12; ++i) T *= 10;
	uint64_t p_max = 1000000;
	if (p_max < 7) p_max = 7;
	if (p_max > uint64_t(-1) / 4) p_max = uint64_t(-1) / 4;

	uint64_t k_min = (argc > 2) ? std::atoll(argv[2]) : 15;
	k_min /= 30; k_min *= 30; k_min += 15;

	uint64_t k_max = (argc > 3) ? std::atoll(argv[3]) : k_min + 30 * uint64_t(143165576);
	// uint64_t k_max = (argc > 3) ? std::atoll(argv[3]) : k_min + 30 * uint64_t(1000000);
	k_max /= 30; k_max *= 30; k_max += 15;
	if (k_max < k_min) k_max = k_min;

	qsieve(p_max, k_min, k_max);

	return EXIT_SUCCESS;
}
