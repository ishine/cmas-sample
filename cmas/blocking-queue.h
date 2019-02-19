#ifndef BLOCKING_QUEUE_H
#define BLOCKING_QUEUE_H

#include <condition_variable>
#include <list>
#include <mutex>
#include <queue>
#include <thread>
#include <memory>
#include <cstddef>
#include <chrono>
#include <iostream>

namespace spicax {
// Special Blocking Queue used for voiceAI voice capture system.
// Each element in the queue is a pointer(i.e. short*, 16 bit), thus users should maintain the memory
// by yourself. all the leftover pointer in the queue will be deleted when calling method Clear()
// or destructing
class BlockingQueue {
 public:
  static constexpr size_t kInfiniteQueueSize = 0;
  // defaultly construct a blocking queue with infinite queue size
  BlockingQueue() : BlockingQueue(kInfiniteQueueSize) {};
  // construct a blocking queue with a size of 'max_size'
  explicit BlockingQueue(const size_t max_size): max_size_(max_size) {};
  // disable copy constructor and assignment operator
  BlockingQueue(const BlockingQueue &) = delete;
  BlockingQueue& operator=(const BlockingQueue &) = delete;
  // destructor, clear all the leftover data
  ~BlockingQueue() {
    Clear();
  };
  // Push an item into the queue, block if it is full
  void Push(short* item) {
    std::unique_lock<std::mutex> lock(mutex_);
    condvar_.wait(lock, [this]() {return NotFull();});
    queue_.push(item);
    condvar_.notify_one();
  };
  // like Push, return nullptr, but currently we do not test this method
  bool PushWithTimeout(short* item, const std::chrono::milliseconds timeout) {
    std::unique_lock<std::mutex> lock(mutex_);
    if (!condvar_.wait_for(lock, timeout, [this]() {return NotFull();})) {
      return false;
    }
    queue_.push(item);
    condvar_.notify_one();
    return true;
  };
  // Pop out an item from the queue, block until an item is available
  short* Pop() {
    std::unique_lock<std::mutex> lock(mutex_);
    condvar_.wait(lock, [this]() { return NotEmpty();});
    short* item = queue_.front();
    queue_.pop();
    return item;
  };
  // like Pop, return nullptr if time out
  short * PopWithTimeOut(const std::chrono::milliseconds timeout) {
    std::unique_lock<std::mutex> lock(mutex_);
    if (!condvar_.wait_for(lock, timeout, [this]() { return NotEmpty();})) {
      return nullptr;
    }
    short* item = queue_.front();
    queue_.pop();
    return item;
  };
  // clear this queue
  void Clear() {
    std::unique_lock<std::mutex> lock(mutex_);
    int qsize = queue_.size();
    for (int i = 0; i < qsize; i++) {
      short *item = queue_.front();
      queue_.pop();
      if (item != nullptr) {delete [] item;}
    }
  };
  int Size() {
    std::unique_lock<std::mutex> lock(mutex_);
    return queue_.size();
  };
 private:
  bool NotFull() const { return max_size_ == kInfiniteQueueSize || queue_.size() < max_size_; }
  bool NotEmpty() const { return !queue_.empty(); }
  const size_t max_size_;
  std::queue<short*, std::list<short*>> queue_;
  mutable std::mutex mutex_;
  std::condition_variable condvar_;
};
}
#endif
