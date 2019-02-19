#ifndef CIRCULAR_BUFFER_H
#define CIRCULAR_BUFFER_H

#include <condition_variable>
#include <memory>
#include <mutex>
#include <cstring>
#include <cstddef>
#include <algorithm>
#include <iostream>
#include "kaldi-common.h"

// make sure (input_size + requried_size) <= size_ to avoid deadlock
namespace spicax {
// template<typename T>
class CircularBuffer {
 public:
  CircularBuffer(): data_(nullptr) {};
  explicit CircularBuffer(size_t size) {Resize(size);};
  CircularBuffer(const CircularBuffer& other) = delete;
  CircularBuffer& operator=(const CircularBuffer & other) = delete;

  ~CircularBuffer() {
    // std::cout << "destructing CircularBuffer\n";
    if (data_ != nullptr) {delete[] data_; data_ = nullptr;}
  }
  // resize circular buffer
  void Resize(size_t size) {
    std::unique_lock<std::mutex> lock(mutex_);
    memory_size_ = size + 1;
    if (data_ != nullptr) {delete[] data_; data_ = nullptr;}
    data_ = new short[memory_size_];
    // data_.reset(new T[memory_size_]);
    head_ = 0;
    tail_ = 0;
    actual_size_ = 0;
    size_ = size;
  };
  // push a block data into data buffer
  void Push(short* buffer, size_t input_size) {
    std::unique_lock<std::mutex> lock(mutex_);
    condvar_.wait(lock, [this, input_size]() {return (actual_size_ + input_size) < memory_size_;});
    size_t part1_size = std::min(head_ + input_size, memory_size_) - head_;
    size_t part2_size = input_size - part1_size;
    std::memcpy(&data_[head_], buffer, sizeof(short) * part1_size);
    std::memcpy(&data_[0], buffer + part1_size, sizeof(short) * part2_size);
    actual_size_ += input_size;
    head_ = (head_ + input_size) % memory_size_;
    condvar_.notify_one();
  };
  //
  void PushWithOverflow(short *buffer, size_t input_size) {
    std::unique_lock<std::mutex> lock(mutex_);
    size_t part1_size = std::min(head_ + input_size, memory_size_) - head_;
    size_t part2_size = input_size - part1_size;
    std::memcpy(&data_[head_], buffer, sizeof(short) * part1_size);
    std::memcpy(&data_[0], buffer + part1_size, sizeof(short) * part2_size);
    actual_size_ += input_size;
    bool overflow = false;
    if (actual_size_ + input_size > memory_size_) {
      overflow = true;
      KALDI_LOG << "input_size=" << input_size << ",actual_size_=" << actual_size_ << ",memory_size_=" << memory_size_ << "\n";
      actual_size_ = size_;
    }
    head_ = (head_ + input_size) % memory_size_;
    if (overflow) {
      tail_ = (head_ + 1) % memory_size_;
    }
    // std::cout << "overflow=" << overflow << ",head_=" << head_ << ",tail_=" << tail_ << "\n";
    condvar_.notify_one();
  }
  // pop out a block of data, buffer is allocated outside
  void Pop(short* buffer, size_t required_size) {
    std::unique_lock<std::mutex> lock(mutex_);
    condvar_.wait(lock, [this, required_size]() {return (actual_size_ >= required_size);});
    size_t part1_size = std::min(tail_ + required_size, memory_size_) - tail_;
    size_t part2_size = required_size - part1_size;
    std::memcpy(buffer, &data_[tail_], sizeof(short) * part1_size);
    std::memcpy(buffer + part1_size, &data_[0], sizeof(short) * part2_size);
    actual_size_ -= required_size;
    tail_ = (tail_ + required_size) % memory_size_;
    condvar_.notify_one();
  };
  //pop out a block of data, with timeout
  int PopWithTimeOut(short *buffer, size_t required_size, const std::chrono::milliseconds timeout) {
    std::unique_lock<std::mutex> lock(mutex_);
    if (!condvar_.wait_for(lock, timeout, [this, required_size]() { return (actual_size_ >= required_size);})) {
      KALDI_LOG << "insufficient data, required_size=" << required_size << ", while actual_size_=" << actual_size_ << "\n";
      size_t part1_size = std::min(tail_ + required_size, memory_size_) - tail_;
      size_t part2_size = required_size - part1_size;
      std::memcpy(buffer, &data_[tail_], sizeof(short) * part1_size);
      std::memcpy(buffer + part1_size, &data_[0], sizeof(short) * part2_size);
      tail_ = head_;
      int actual_size = actual_size_;
      actual_size_ = 0;
      return actual_size;
    }
    size_t part1_size = std::min(tail_ + required_size, memory_size_) - tail_;
    size_t part2_size = required_size - part1_size;
    std::memcpy(buffer, &data_[tail_], sizeof(short) * part1_size);
    std::memcpy(buffer + part1_size, &data_[0], sizeof(short) * part2_size);
    actual_size_ -= required_size;
    tail_ = (tail_ + required_size) % memory_size_;
    return required_size;
    condvar_.notify_one();
  }
  // Set buffer to empty but not resize
  void Clear() {
    std::unique_lock<std::mutex> lock(mutex_);
    tail_ = head_;
  };
  size_t Size() {
    std::unique_lock<std::mutex> lock(mutex_);
    return size_;
  }

 private:
  std::mutex mutex_;
  std::condition_variable condvar_;
  // std::unique_ptr<T[]> data_;
  short *data_;
  // Full: head_ = (tail_ + 1) % memory_size_; Empty: head_ = tail_;
  size_t head_;//
  size_t tail_;//
  size_t size_;//
  size_t actual_size_;//
  size_t memory_size_;//
};

// template class CircularBuffer<short>;
}
#endif