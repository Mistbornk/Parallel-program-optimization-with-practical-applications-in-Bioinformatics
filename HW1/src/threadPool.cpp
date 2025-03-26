#include "threadPool.hpp"
#include <iostream>
#include <sstream>
std::vector<std::string> thread_logs; 

ThreadPool::ThreadPool(size_t num_threads) : stop(false) {
	if (num_threads < 5) {
		num_threads = 5;
		std::cout << "Must at least 5 thread ! Assign thread number equal to 5." << std::endl;
	}
	auto work = [&](void) {
		auto thread_id = std::this_thread::get_id();
		auto start_time = std::chrono::steady_clock::now();
		std::chrono::duration<double> total_runtime(0);
		
		while (true) {
			std::function<void()> job;
			// scope of mutex lock
			{
				std::unique_lock<std::mutex> lock(mtx); 
				cv.wait(lock, [this] {return stop || !jobs.empty();}); // while break, unlock mutex;
				// if all jobs complete and no more job 
				if (stop && jobs.empty()) break;
				// select a job from queue
				job = std::move(jobs.front()); 
				jobs.pop(); 
			}
			// do the job
			auto job_start_time = std::chrono::steady_clock::now();
			job();
			auto job_end_time = std::chrono::steady_clock::now();
			total_runtime += job_end_time - job_start_time;
		}
		auto end_time = std::chrono::steady_clock::now();
		std::chrono::duration<double> total_lifetime = end_time - start_time;
		// store info
        {
            std::lock_guard<std::mutex> lock(mtx);
            std::ostringstream log;
            log << "\n------------------\n";
            log << "Thread id: " << thread_id << "\n"
                << "total runtime:  " << total_runtime.count() << "s\n"
                << "total lifetime: " << total_lifetime.count() << "s\n";
            thread_logs.push_back(log.str());
        }
	};
	for (size_t t=0; t<num_threads; t++) {
		threads.emplace_back(work);
	}
}

void ThreadPool::enqueueJobs(std::function<void()> job) {
	// scope of mutex lock
	{
		std::lock_guard<std::mutex> lock(mtx);
		jobs.push(job);
	}
	cv.notify_one();
}

ThreadPool::~ThreadPool() {
	// scope of mutex lock
	{
		std::lock_guard<std::mutex> lock(mtx);
		stop = true;
	}
	cv.notify_all();
	// join thread
	for (std::thread &thread : threads) {
		if (thread.joinable()) thread.join();
	}
	std::cout << "\n\n[ThreadPool Summary]\n";
	for (const auto &log : thread_logs) {
		std::cout << log;
	}
}