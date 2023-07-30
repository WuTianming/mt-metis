/**
 * prefetch functions and structures, local to this file
*/

#define GNU_SOURCE
#include <stdio.h>
#include <fcntl.h>
#include <liburing.h>

const static uint8_t CACHE_INVALID = 0;
const static uint8_t CACHE_SHARED = 1;
const static uint8_t CACHE_FETCHING = 2;
const static uint8_t CACHE_TODISCARD = 3;
static struct pagecacher {
  int64_t pagesize;
  int64_t pagecnt;
  uint64_t capacity;
  uint64_t resident;
  int64_t **cached_pages;

  int fd;
  struct io_uring ring;
  struct io_uring_cqe *cqe; // completion queue entry

  struct pageread {
    int64_t page_id;
    struct iovec iov;
  };
  uint8_t *cached_status;
  int64_t m;
};

static void PageCacher(struct pagecacher *this,
      uint64_t capacity, int64_t m, int64_t pagesize, int fd) {
  this->capacity = capacity;  // unused for now
  this->resident = 0;
  this->pagesize = pagesize;
  this->m = m;
  this->pagecnt = (m - 1) / pagesize + 1;

  this->cached_pages = (int64_t **)malloc(sizeof(int64_t *) * this->pagecnt);
  this->cached_status = (uint8_t *)malloc(sizeof(uint8_t) * this->pagecnt);

  memset(this->cached_pages, 0, sizeof(void *) * this->pagecnt);
  memset(this->cached_status, CACHE_INVALID, sizeof(uint8_t) * this->pagecnt);

  this->fd = fd;

  io_uring_queue_init(256, &(this->ring), 0);
}

static void PageCacher_Del(struct pagecacher *this) {
  io_uring_queue_exit(&this->ring);
  for (int i = 0; i < this->pagecnt; i++)
    if (this->cached_status[i] != CACHE_INVALID)
      free(this->cached_pages[i]);
  free(this->cached_pages);
  close(this->fd);
}

#include <stdbool.h>
static bool prefetch_page(struct pagecacher *this, int64_t p, bool force) {
  // send an request using linux liburing
  // store data in cached_pages
  int64_t pagecnt = this->pagecnt;
  int64_t pagesize = this->pagesize;
  if (p >= pagecnt)
    return true;
  if (!force && this->resident >= this->capacity)
    return false; // only returns false when capacity is reached
  if (this->cached_status[p] == CACHE_FETCHING ||
      this->cached_status[p] == CACHE_SHARED)
    return true;
  if (this->cached_status[p] == CACHE_TODISCARD) {
    // cancel discard and proceed with fetch
    this->cached_status[p] = CACHE_FETCHING;
    this->resident++;
    return true;
  }
  this->cached_status[p] = CACHE_FETCHING;
  this->resident++;

  size_t len = pagesize * sizeof(int64_t);
  int64_t page_offset = p * len;

  struct pageread *pr = malloc(sizeof(struct pageread));
  pr->page_id = p;
  struct iovec *iov = &pr->iov;
  iov->iov_len = len;
  int ret = posix_memalign((void **)&iov->iov_base, 4096, len);
  if (iov->iov_base == NULL) {
    perror("posix_memalign");
    exit(1);
  }

  // this->cached_pages[p] = (int64_t *)iov->iov_base;

  struct io_uring_sqe *sqe = io_uring_get_sqe(&this->ring);
  io_uring_prep_readv(sqe, this->fd, iov, 1, page_offset);
  io_uring_sqe_set_data(sqe, pr);
  ret = io_uring_submit(&this->ring);

  if (ret < 1) {
    perror("io_uring_submit");
    exit(1);
  }
  return true;
}

static void process_uring(struct pagecacher *this) {
  // process the completion queue
  int ret = io_uring_wait_cqe(&this->ring, &this->cqe);
  if (ret < 0) {
    perror("io_uring_wait_cqe");
    exit(1);
  }
  if (this->cqe->res < 0) {
    fprintf(stderr, "Async readv failed.\n");
    exit(1);
  }

  struct pageread *pr = (struct pageread *)io_uring_cqe_get_data(this->cqe);
  int64_t page_id = pr->page_id;
  int64_t expected = this->pagesize * sizeof(int64_t);

  if (page_id == this->pagecnt - 1)
    expected = (this->m % this->pagesize) * sizeof(int64_t);
  if (this->cqe->res != expected) {
    fprintf(stderr,
            "On page %ld, async readv reads %d bytes, expected %ld bytes.\n",
            page_id, this->cqe->res, expected);
    exit(1);
  }

  if (this->cached_status[page_id] == CACHE_TODISCARD) {
    this->cached_status[page_id] = CACHE_INVALID;
    free(pr->iov.iov_base); // discard read data
  } else if (this->cached_status[page_id] == CACHE_FETCHING) {
    this->cached_status[page_id] = CACHE_SHARED;
    this->cached_pages[page_id] = (int64_t *)pr->iov.iov_base;
  }

  io_uring_cqe_seen(&this->ring, this->cqe);
  free(pr);
}

static void fetch_page(struct pagecacher *this, int64_t p) {
  if (p >= this->pagecnt)
    return;

  prefetch_page(this, p, true);

  while (this->cached_status[p] != CACHE_SHARED) {
    process_uring(this);
  }
}

static void discard_page(struct pagecacher *this, int64_t p) {
  if (p >= this->pagecnt)
    return;
  switch (this->cached_status[p]) {
  case CACHE_INVALID:
    break;
  case CACHE_SHARED:
    this->cached_status[p] = CACHE_INVALID;
    free(this->cached_pages[p]);
    this->resident--;
    break;
  case CACHE_FETCHING:
    this->cached_status[p] = CACHE_TODISCARD;
    this->resident--;
    // discard query will be processed when fetch is done
    break;
  case CACHE_TODISCARD:
    break;
  }
  this->cached_pages[p] = NULL;
}