#include <iostream>
#include "threadPool.hpp"
#include <atomic>
#include <mutex>
#include <condition_variable>

using namespace std;

std::mutex cout_mutex;
std::mutex cv_mutex;
std::condition_variable cv;
std::atomic<int> pending_jobs(496);  // avoid race condition

void print_1() {
    int num = rand() % 100;
    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        cout << ((num % 2) ? '1' : '0');
    }
    if (--pending_jobs == 0) {  // notice condictional varialble when pending_jobs==0 
        std::lock_guard<std::mutex> lock(cv_mutex);
		cout << "\n\nprint_2: " << endl;
        cv.notify_all();
    }
}

struct Print2 {
    void operator()() {
        std::unique_lock<std::mutex> lock1(cv_mutex);
        cv.wait(lock1, [] { return pending_jobs == 0; });

		std::lock_guard<std::mutex> lock2(cout_mutex); 
		std::cout << "2";
    }
};

int main() {
    ThreadPool pool(5);
	cout << "print_1: " << endl;
    for (int i = 0; i < 496; ++i) {
        pool.enqueueJobs(print_1);
    }
    for (int i = 0; i < 4; ++i) {
        pool.enqueueJobs(Print2());
    }
    return 0;
}
