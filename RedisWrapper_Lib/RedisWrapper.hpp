#ifndef REDISWRAPPER_HPP
#define REDISWRAPPER_HPP

#include <mutex>

#include "hiredis.h"

class RedisWrapper
{
public:
    static RedisWrapper& getContext(const char * redisHost = "localhost");

    redisReply * redisCommand(const char *format, ...);

    static const unsigned int MAX_HOST_LENGTH = 64;

private:
    RedisWrapper();
    RedisWrapper(const char * redisHost);
    ~RedisWrapper();

    redisContext * rContext;

    bool hostSet;
    std::mutex rMutex;

    void setHost(const char * redisHost);
    RedisWrapper(RedisWrapper const&) = delete;
    void operator=(RedisWrapper const&) = delete;
};

#endif
