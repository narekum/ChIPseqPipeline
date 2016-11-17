# From http://code.activestate.com/recipes/498245-lru-and-lfu-cache-decorators/
# McBryan: Added Destructor
# McBryan: Remove old cache entries then add new one, fixes issue where cache entry is added then immediately removed as the least frequent

import collections
import functools
from itertools import ifilterfalse
from heapq import nsmallest
from operator import itemgetter

class Counter(dict):
    'Mapping where default values are zero'
    def __missing__(self, key):
        return 0

def lru_cache(maxsize=100, destructor = None):
    '''Least-recently-used cache decorator.

    Arguments to the cached function must be hashable.
    Cache performance statistics stored in f.hits and f.misses.
    Clear the cache with f.clear().
    http://en.wikipedia.org/wiki/Cache_algorithms#Least_Recently_Used

    '''
    maxqueue = maxsize * 10
    def decorating_function(user_function,
            len=len, iter=iter, tuple=tuple, sorted=sorted, KeyError=KeyError):
        cache = {}                  # mapping of args to results
        queue = collections.deque() # order that keys have been used
        refcount = Counter()        # times each key is in the queue
        sentinel = object()         # marker for looping around the queue
        kwd_mark = object()         # separate positional and keyword args

        # lookup optimizations (ugly but fast)
        queue_append, queue_popleft = queue.append, queue.popleft
        queue_appendleft, queue_pop = queue.appendleft, queue.pop

        @functools.wraps(user_function)
        def wrapper(*args, **kwds):
            # cache key records both positional and keyword args
            key = args
            if kwds:
                key += (kwd_mark,) + tuple(sorted(kwds.items()))

            # get cache entry or compute if not found
            try:
                result = cache[key]
                wrapper.hits += 1
            except KeyError:
                # purge least recently used cache entry
                if len(cache) >= maxsize:
                    dkey = queue_popleft()
                    refcount[dkey] -= 1
                    # remove any duplicate entries in the queue
                    while refcount[dkey]:
                        dkey = queue_popleft()
                        refcount[dkey] -= 1
                    if destructor != None:
                        destructor(cache[dkey])
                    del cache[dkey], refcount[dkey]
                    
                # record recent use of this key
                queue_append(key)
                refcount[key] += 1
                    
                result = user_function(*args, **kwds)
                cache[key] = result
                wrapper.misses += 1

            # periodically compact the queue by eliminating duplicate keys
            # while preserving order of most recent access
            if len(queue) > maxqueue:
                refcount.clear()
                queue_appendleft(sentinel)
                for key in ifilterfalse(refcount.__contains__,
                                        iter(queue_pop, sentinel)):
                    queue_appendleft(key)
                    refcount[key] = 1


            return result

        def clear():
            if destructor != None:
                for key in cache:
                    destructor(cache[key])
            cache.clear()
            queue.clear()
            refcount.clear()
            wrapper.hits = wrapper.misses = 0

        wrapper.hits = wrapper.misses = 0
        wrapper.clear = clear
        return wrapper
    return decorating_function




def lfu_cache(maxsize=100, destructor = None):
    '''Least-frequenty-used cache decorator.

    Arguments to the cached function must be hashable.
    Cache performance statistics stored in f.hits and f.misses.
    Clear the cache with f.clear().
    http://en.wikipedia.org/wiki/Least_Frequently_Used

    '''
    def decorating_function(user_function):
        cache = {}                      # mapping of args to results
        use_count = Counter()           # times each key has been accessed
        kwd_mark = object()             # separate positional and keyword args

        @functools.wraps(user_function)
        def wrapper(*args, **kwds):
            key = args
            if kwds:
                key += (kwd_mark,) + tuple(sorted(kwds.items()))
            
            # get cache entry or compute if not found
            try:
                result = cache[key]
                wrapper.hits += 1
            except KeyError:

                # purge least frequently used cache entries
                if len(cache) >= maxsize:
                    
                    # inner join for (key,count) for only keys which are currently in the cache
                    # this means we can keep a forever log of keys we've found and remove based on long term history
                    current_key_counts = [(r1,use_count[r2]) for r1 in cache.keys() for r2 in use_count.keys() if r1==r2]
                    
                    for dkey, _ in nsmallest(maxsize // 10,
                                            current_key_counts,
                                            key=itemgetter(1)):
                        if destructor != None:
                            destructor(cache[dkey])
                        del cache[dkey]
                
                use_count[key] += 1       
                result = user_function(*args, **kwds)
                cache[key] = result
                wrapper.misses += 1

            return result

        def clear():
            if destructor != None:
                for key in cache:
                    destructor(cache[key])
            cache.clear()
            use_count.clear()
            wrapper.hits = wrapper.misses = 0

        wrapper.hits = wrapper.misses = 0
        wrapper.clear = clear
        return wrapper
    return decorating_function


if __name__ == '__main__':

    @lru_cache(maxsize=20)
    def f(x, y):
        return 3*x+y

    domain = range(5)
    from random import choice
    for i in range(1000):
        r = f(choice(domain), choice(domain))

    print(f.hits, f.misses)

    @lfu_cache(maxsize=20)
    def f(x, y):
        return 3*x+y

    domain = range(5)
    from random import choice
    for i in range(1000):
        r = f(choice(domain), choice(domain))

    print(f.hits, f.misses)