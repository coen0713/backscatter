#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace bsar {

/// 0 means "use every hardware thread".
inline unsigned resolve_threads(unsigned requested) {
  if (requested != 0) {
    return requested;
  }
  const unsigned hw = std::thread::hardware_concurrency();
  return hw == 0 ? 1 : hw;
}

/// Run `fn(i)` for every i in [begin, end) on up to `threads` threads, handing
/// out chunks of `grain` indices dynamically. Callers keep results
/// deterministic by making each index write only its own output; scheduling
/// order then has no effect on the result.
template <typename Fn>
void parallel_for(std::size_t begin, std::size_t end, unsigned threads, Fn&& fn,
                  std::size_t grain = 1) {
  if (end <= begin) {
    return;
  }
  grain = std::max<std::size_t>(grain, 1);
  const std::size_t chunks = (end - begin + grain - 1) / grain;
  const auto workers =
      static_cast<unsigned>(std::min<std::size_t>(resolve_threads(threads), chunks));
  if (workers <= 1) {
    for (std::size_t i = begin; i < end; ++i) {
      fn(i);
    }
    return;
  }

  std::atomic<std::size_t> next{begin};
  std::exception_ptr error;
  std::mutex error_mutex;
  auto work = [&] {
    try {
      for (;;) {
        const std::size_t start = next.fetch_add(grain, std::memory_order_relaxed);
        if (start >= end) {
          return;
        }
        const std::size_t stop = std::min(end, start + grain);
        for (std::size_t i = start; i < stop; ++i) {
          fn(i);
        }
      }
    } catch (...) {
      const std::lock_guard<std::mutex> lock(error_mutex);
      if (!error) {
        error = std::current_exception();
      }
      next.store(end, std::memory_order_relaxed);
    }
  };

  std::vector<std::thread> pool;
  pool.reserve(workers - 1);
  for (unsigned t = 1; t < workers; ++t) {
    pool.emplace_back(work);
  }
  work();
  for (auto& th : pool) {
    th.join();
  }
  if (error) {
    std::rethrow_exception(error);
  }
}

}  // namespace bsar
