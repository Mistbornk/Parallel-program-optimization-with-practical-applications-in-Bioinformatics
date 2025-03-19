#ifndef THREAD_POOL_H
#define THREAD_POOL_H
#include <queue>
#include <vector>
#include <thread>
#include <functional>
#include <condition_variable>
#include <chrono>
#include <atomic>

class ThreadPool {
public:
	ThreadPool(size_t num_threads = 5);
	~ThreadPool();
	void enqueueJobs(std::function<void()> job);

private:
	std::vector<std::thread> threads;
	std::queue<std::function<void()>> jobs;
	std::mutex mtx;
	std::condition_variable cv;
	bool stop;
};

#endif