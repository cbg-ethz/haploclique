/* Copyright 2012 Tobias Marschall
 * 
 * This file is part of HaploClique.
 * 
 * HaploClique is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HaploClique is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HaploClique.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <iostream>
#include <memory>
#include <deque>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <queue>
#include <limits>

#include <pthread.h>
#include <semaphore.h>

template <typename WorkPackageType, typename OutputWriterType>
class ThreadPool {
private:
	typedef struct job_t {
		std::vector<WorkPackageType*> work_packages;
		long long id;
		job_t(long long id) : id(id) {}
	} job_t;
	typedef struct job_comparator_t {
		bool operator()(const job_t* e1, const job_t* e2) { return e1->id > e2->id; }
	} job_comparator_t;
	typedef std::priority_queue<job_t*, std::vector<job_t*>, job_comparator_t> output_queue_t;

	int worker_thread_count;
	int job_size;
	OutputWriterType& output_writer;
	std::deque<job_t*> queue;
	job_t* new_job;
	
	pthread_t* threads;
	sem_t queue_add_semaphore;
	sem_t queue_remove_semaphore;
	pthread_mutex_t queue_mutex;
	long long queued_job_count;
	long long finished_job_count;
	
	output_queue_t output_queue;
	pthread_mutex_t output_mutex;
	pthread_t output_thread;
	
	static void *thread_main(void *p) {
		assert(p != 0);
		typedef ThreadPool<WorkPackageType,OutputWriterType> thread_pool_t;
		thread_pool_t* parent = (thread_pool_t*)p;
		while (true) {
			sem_wait(&parent->queue_remove_semaphore);
			pthread_mutex_lock(&parent->queue_mutex);
			assert(parent->queue.size() > 0);
			job_t* e = parent->queue[0];
			parent->queue.pop_front();
			sem_post(&parent->queue_add_semaphore);
			pthread_mutex_unlock(&parent->queue_mutex);
			assert(e != 0);
			// was termination of thread requested?
			if (e->id == -1) {
				delete e;
				e = 0;
			} else {
				// do the actual work
				for (size_t i=0; i<e->work_packages.size(); ++i) {
					assert(e->work_packages[i] != 0);
					e->work_packages[i]->run();
				}
			}
			// write output
			pthread_mutex_lock(&parent->output_mutex);
			if (e != 0) {
				parent->output_queue.push(e);
			}
			while (!parent->output_queue.empty()) {
				job_t* job = parent->output_queue.top();
				assert(job != 0);
				if (job->id == parent->finished_job_count) {
					for (size_t j=0; j<job->work_packages.size(); ++j) {
						assert(job->work_packages[j] != 0);
						parent->output_writer.write(std::auto_ptr<WorkPackageType>(job->work_packages[j]));
					}
					parent->output_queue.pop();
					delete job;
					parent->finished_job_count += 1;
				} else {
					break;
				}
			}
			pthread_mutex_unlock(&parent->output_mutex);
			if (e == 0) break;
		}
		return NULL;
	}

	void queue_job(job_t* job) {
		sem_wait(&queue_add_semaphore);
		pthread_mutex_lock(&queue_mutex);
		assert(job != 0);
		assert(job->id == queued_job_count);
		queue.push_back(job);
		queued_job_count += 1;
		pthread_mutex_unlock(&queue_mutex);
		sem_post(&queue_remove_semaphore);
	}
	
public:
	/** If worker_thread_count equals 0, no additional threads will be spawned and all operations are performed in the invoking thread. */
	ThreadPool(int worker_thread_count, int job_size, int queue_size, OutputWriterType& output_writer) : worker_thread_count(worker_thread_count), job_size(job_size), output_writer(output_writer) {
		assert(worker_thread_count >= 0);
		assert(job_size > 0);
		queued_job_count = 0;
		finished_job_count = 0;
		if (worker_thread_count > 0) {
			new_job = new job_t(0);
			if (sem_init(&queue_add_semaphore, 0, queue_size)) {
				throw std::runtime_error("Could not create semaphore");
			}
			if (sem_init(&queue_remove_semaphore, 0, 0)) {
				throw std::runtime_error("Could not create semaphore");
			}
			if (pthread_mutex_init(&queue_mutex, NULL)) {
				throw std::runtime_error("Could not create mutex");
			}
			if (pthread_mutex_init(&output_mutex, NULL)) {
				throw std::runtime_error("Could not create mutex");
			}
			threads = new pthread_t[worker_thread_count];
			for (size_t i=0; i<worker_thread_count; ++i) {
				if (pthread_create(threads+i, NULL, &thread_main, this)) {
					throw std::runtime_error("Could not create thread");
				}
			}
		}
	}
	
	~ThreadPool() {
		if (worker_thread_count > 0) {
			pthread_mutex_lock(&queue_mutex);
			if (new_job->work_packages.size() > 0) {
				queue.push_back(new_job);
				sem_post(&queue_remove_semaphore);
				new_job = 0;
			}
			for (int i=0; i<worker_thread_count; ++i) {
				queue.push_back(new job_t(-1));
				sem_post(&queue_remove_semaphore);
			}
			pthread_mutex_unlock(&queue_mutex);
			for (int i=0; i<worker_thread_count; ++i) {
				if (pthread_join(threads[i], NULL)) {
					throw std::runtime_error("Error during phtread_join");
				}
			}
			delete [] threads;
			sem_destroy(&queue_add_semaphore);
			sem_destroy(&queue_remove_semaphore);
			pthread_mutex_destroy(&queue_mutex);
			pthread_mutex_destroy(&output_mutex);
			if (new_job != 0) delete new_job;
		}
	}
	
	void addTask(std::auto_ptr<WorkPackageType> work_package) {
		assert(work_package.get() != 0);
		if (worker_thread_count == 0) {
			work_package->run();
			output_writer.write(work_package);
		} else {
			assert(new_job != 0);
			new_job->work_packages.push_back(work_package.release());
			if (new_job->work_packages.size() == job_size) {
				queue_job(new_job);
				new_job = new job_t(queued_job_count);
			}
		}
	}
	
};

#endif /* THREADPOOL_H_ */
